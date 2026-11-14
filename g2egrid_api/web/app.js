"use strict";
/* g2egrid -- interface web.
 *
 * Un seul objet d'etat, E, qui porte LES TABLES du reseau : c'est la meme
 * description que les notebooks et le MATLAB, pas une representation
 * intermediaire. Tout ce qu'on dessine se traduit immediatement en lignes de
 * table (onglet Tables pour s'en convaincre), et toutes les routes de calcul
 * recoivent ces tables en corps de requete.
 *
 * E.version est incremente a chaque edition. Chaque resultat retient la
 * version pour laquelle il a ete calcule, et toute reponse asynchrone verifie
 * cette version avant de s'afficher : editer pendant un calcul annule
 * silencieusement l'affichage du resultat perime, au lieu de montrer une carte
 * des modes qui ne correspond plus au reseau.
 */

const E = {
  tables: null,
  positions: {},
  defauts: null,
  version: 0,
  outil: "select",
  sel: null, // {table, index}
  attente: null, // premier bus d'une branche en cours
  drag: null,
  // Un reseau a zones se tient sous sa forme STRUCTUREE -- transport, zones,
  // instances -- et le constructeur en edite une VUE a la fois : le
  // transport, ou une zone. E.tables pointe alors sur les tables de la vue,
  // si bien que tous les outils du constructeur marchent sans le savoir.
  // Les calculs, eux, prennent E.plat : le reseau aplati par le serveur.
  reseau: null,
  vue: null,
  plat: null,
  erreurAplat: "",
  // indice DANS LE SPECTRE du mode ouvert sous la carte
  mode: null,
  // le scenario en cours d'edition (onglet Temporel)
  scenario: [],
  // hypotheses de simplification (onglet Simplifications)
  simp: {},
  simpAppliquer: false,
  simpInv: null,
  simpInvVersion: -1,
  res: { lf: null, modal: null, signaux: null, forme: null, balayage: null, validation: null },
};

// =====================================================================
// utilitaires
// =====================================================================
const $ = (s) => document.querySelector(s);
const $$ = (s) => Array.from(document.querySelectorAll(s));

function el(nom, attrs = {}, ...enfants) {
  const svg = ["svg", "g", "circle", "line", "path", "text", "rect", "polyline", "polygon"];
  const n = svg.includes(nom)
    ? document.createElementNS("http://www.w3.org/2000/svg", nom)
    : document.createElement(nom);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") n.setAttribute("class", v);
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined) n.setAttribute(k, v);
  }
  for (const c of enfants) if (c !== null && c !== undefined) n.append(c);
  return n;
}

const fmt = (v, n = 4) =>
  v === null || v === undefined || Number.isNaN(v)
    ? "--"
    : Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0)
      ? v.toExponential(2)
      : (+v).toFixed(n).replace(/\.?0+$/, "");

function etat(msg, erreur = false) {
  const n = $("#etat");
  n.textContent = msg || "";
  n.style.color = erreur ? "var(--erreur)" : "var(--doux)";
}

function modifie() {
  E.version += 1;
  E.res = { lf: null, modal: null, signaux: null, forme: null, balayage: null, validation: null, impact: null };
  E.disj = null;
  // le mode ouvert etait un indice dans l'ANCIEN spectre : le garder
  // montrerait une forme modale qui ne decrit plus le reseau affiche
  E.mode = null;
  const rf = $("#rangee-forme");
  if (rf) rf.hidden = true;
  const rs = $("#rangee-mode-schema");
  if (rs) rs.hidden = true;
  arreterAnimationMode();
  const rr = $("#rangee-rejeu");
  if (rr) rr.hidden = true;
  arreterRejeu();
  E._rejeuT = null;
  if (E.reseau) aplatirBientot();
  dessiner();
  rendreTables();
  rendrePoint();
  validerBientot();
  sauverBrouillon();
  memoriser(true);
}

// =====================================================================
// annuler / retablir
// =====================================================================
// L'historique garde des INSTANTANES complets du reseau -- ce que
// l'enregistrement ecrirait --, pas des operations inverses : chaque outil du
// constructeur reste libre d'editer les tables comme il l'entend, et
// l'annulation ne depend d'aucun d'eux. Un instantane est pris a chaque
// modifie() ; des editions rapprochees (une saisie au clavier, les crans d'un
// curseur) se fondent en une seule entree. Revenir a un instantane repasse
// par poser() / entrerZones(), les memes chemins qu'une ouverture de fichier.
const H = { pile: [], i: -1, t: 0, max: 100, gel: false };

function instantane() {
  return JSON.stringify(E.reseau ? { zones: corpsStructure(), vue: E.vue } : { plat: corpsFichier() });
}

function memoriser(fusion = false) {
  if (H.gel || !E.tables) return;
  const s = instantane();
  if (H.pile[H.i] === s) return;
  const maintenant = Date.now();
  H.pile = H.pile.slice(0, H.i + 1);
  if (fusion && H.i > 0 && maintenant - H.t < 600) H.pile[H.i] = s;
  else {
    H.pile.push(s);
    if (H.pile.length > H.max) H.pile.shift();
  }
  H.i = H.pile.length - 1;
  H.t = maintenant;
  majHistorique();
}

function majHistorique() {
  const a = $("#btn-annuler"), r = $("#btn-retablir");
  if (a) a.disabled = H.i <= 0;
  if (r) r.disabled = H.i >= H.pile.length - 1;
}

function revenir(pas) {
  const j = H.i + pas;
  if (j < 0 || j >= H.pile.length) return;
  const s = JSON.parse(H.pile[j]);
  H.gel = true;
  try {
    if (s.zones) {
      entrerZones({ ...s.zones, tables: null, positions: null });
      E.plat = null;
      aplatirBientot();
      if (s.vue && s.vue !== "transport" && E.reseau.zones[s.vue]) ouvrirVue(s.vue);
      dessiner();
      rendreTables();
    } else {
      sortirZones();
      poser(s.plat.tables, s.plat.positions, s.plat.nom);
    }
  } finally {
    H.gel = false;
  }
  H.i = j;
  H.t = 0; // l'edition suivante ouvre une nouvelle entree
  majHistorique();
  // l'onglet Blocs montre le schema de l'instantane retrouve (sous-systemes compris)
  if (typeof rendreBlocs === "function" && $("#onglet-blocs").classList.contains("actif")) rendreBlocs();
  if (typeof rafraichirEquations === "function") rafraichirEquations();
  etat(pas < 0 ? `modification annulee (${H.i} en arriere possibles)` : "modification retablie");
}

// =====================================================================
// brouillon
// =====================================================================
// Un onglet ferme, un navigateur qui plante, un F5 de trop : sans brouillon,
// le reseau dessine depuis une heure est perdu. Chaque edition le range donc
// dans le stockage LOCAL du navigateur (rien ne part au serveur), et la page
// le rouvre au demarrage. Ce n'est pas un enregistrement : le fichier
// .g2egrid.json et le catalogue restent les seules sauvegardes durables -- le
// stockage du navigateur peut etre vide, plein ou interdit (navigation
// privee), et tout le reste doit marcher sans lui.
const CLE_BROUILLON = "g2egrid.brouillon";

function sauverBrouillon() {
  if (E._restauration || !E.tables) return;
  clearTimeout(E._tBrouillon);
  E._tBrouillon = setTimeout(() => {
    try {
      const b = E.reseau
        ? { format: 1, date: Date.now(), reseau: corpsStructure(), vue: E.vue }
        : { format: 1, date: Date.now(), ...corpsFichier() };
      localStorage.setItem(CLE_BROUILLON, JSON.stringify(b));
    } catch {
      // stockage plein ou interdit : le brouillon est un confort, pas un devoir
    }
  }, 600);
}

function lireBrouillon() {
  try {
    const b = JSON.parse(localStorage.getItem(CLE_BROUILLON) || "null");
    return b && b.format === 1 && (b.tables || b.reseau) ? b : null;
  } catch {
    return null;
  }
}

/** Rouvrir le dernier brouillon ; faux s'il n'y en a pas. */
function restaurerBrouillon() {
  const b = lireBrouillon();
  if (!b) return false;
  E._restauration = true;
  try {
    if (b.reseau) {
      entrerZones({ ...b.reseau, tables: null, positions: null });
      E.plat = null;
      aplatirBientot();
      if (b.vue && b.vue !== "transport" && E.reseau.zones[b.vue]) ouvrirVue(b.vue);
    } else {
      sortirZones();
      poser(b.tables, b.positions, b.nom);
    }
  } catch (err) {
    E._restauration = false;
    return false;
  }
  E._restauration = false;
  const quand = new Date(b.date).toLocaleString(I18N.locale, { dateStyle: "short", timeStyle: "short" });
  const nom = (b.reseau && b.reseau.nom) || b.nom || "sans nom";
  etat(`brouillon du ${quand} rouvert (${nom}) -- Enregistrer pour le garder. `);
  $("#etat").append(el("a", {
    href: "#",
    onclick: (ev) => {
      ev.preventDefault();
      effacerBrouillon();
      chargerCatalogue($("#catalogue").value);
    },
  }, "Repartir du catalogue"));
  return true;
}

function effacerBrouillon() {
  try {
    localStorage.removeItem(CLE_BROUILLON);
  } catch {
    // rien a effacer
  }
}

// =====================================================================
// zones
// =====================================================================
function aplatirBientot() {
  E.plat = null; // perime : les calculs attendront le neuf
  clearTimeout(E._tAplat);
  E._tAplat = setTimeout(() => aplatirMaintenant().catch(() => {}), 300);
}

/** Le reseau structure, sous la forme de lecture que le serveur accepte. */
function corpsStructure() {
  const r = E.reseau;
  return { nom: $("#nom-reseau").value || r.nom, f_zoom: r.f_zoom, transport: r.transport, zones: r.zones, instances: r.instances };
}

async function aplatirMaintenant() {
  if (!E.reseau) return;
  clearTimeout(E._tAplat);
  const v = E.version;
  try {
    const d = await api("/api/zones/aplatir", { reseau: corpsStructure() });
    if (v !== E.version) return; // une edition est arrivee entre-temps
    E.plat = { tables: d.tables, positions: d.positions || {} };
    E.erreurAplat = "";
  } catch (err) {
    if (v !== E.version) return;
    E.plat = null;
    E.erreurAplat = err.message;
    throw err;
  }
}

/** Attendre un reseau aplati a jour ; echouer avec le message du serveur si
 *  la structure est incoherente, plutot que de calculer sur la vue seule. */
async function pret() {
  if (!E.reseau || E.plat) return;
  await aplatirMaintenant();
  if (!E.plat) throw new Error(E.erreurAplat || "reseau a zones incoherent");
}

/** Entrer en mode zones a partir d'une reponse de lecture. */
function entrerZones(d) {
  E.reseau = { nom: d.nom, f_zoom: d.f_zoom, transport: d.transport, zones: d.zones, instances: d.instances };
  E.plat = { tables: d.tables, positions: d.positions };
  E.erreurAplat = "";
  $("#nom-reseau").value = d.nom || "";
  adopterSimplifications(d.tables);
  E.version += 1;
  E.res = { lf: null, modal: null, signaux: null, forme: null, balayage: null, validation: null };
  E.disj = null;
  ouvrirVue("transport");
  sauverBrouillon();
  memoriser();
  H.t = 0;
}

function sortirZones() {
  E.reseau = null;
  E.vue = null;
  E.plat = null;
  rendreFilAriane();
}

/** Ouvrir une vue : le transport, ou une zone. Changer de vue ne modifie pas
 *  le reseau -- les resultats deja calcules restent valables. */
function ouvrirVue(v) {
  const cible = v === "transport" ? E.reseau.transport : E.reseau.zones[v];
  if (!cible) return;
  E.vue = v;
  E.tables = cible.tables;
  if (!cible.positions) cible.positions = {};
  E.positions = cible.positions;
  E.sel = null;
  E.attente = null;
  rendreFilAriane();
  dessiner();
  rendreTables();
  validerBientot();
}

function rendreFilAriane() {
  const n = $("#fil-ariane");
  if (!n) return;
  n.textContent = "";
  n.hidden = !E.reseau;
  if (!E.reseau) return;
  const r = E.reseau;
  n.append(el("b", {}, r.nom || "reseau"), " \u203a ");
  if (E.vue === "transport") {
    n.append(el("span", {}, "transport"));
    const nz = Object.keys(r.zones).length;
    n.append(el("span", { class: "doux" },
      `  --  ${nz} zone${nz > 1 ? "s" : ""}, ${r.instances.length} instance${r.instances.length > 1 ? "s" : ""}. ` +
      "Double-cliquer une zone pour l'ouvrir."));
    return;
  }
  const lien = el("a", { href: "#" }, "transport");
  lien.addEventListener("click", (ev) => {
    ev.preventDefault();
    ouvrirVue("transport");
  });
  const z = r.zones[E.vue];
  const inst = r.instances.filter((i) => i.zone === E.vue);
  n.append(lien, " \u203a ", el("span", {}, `zone ${E.vue}`));
  n.append(el("span", { class: "doux" },
    `  --  ${inst.length} instance${inst.length > 1 ? "s" : ""} (bus amont ` +
    inst.map((i) => `${i.bus_amont}${i.K > 1 ? " x" + i.K : ""}`).join(", ") +
    `). Raccordee au transport par son bus ${z.bus_pcc} (PCC). ` +
    "Une modification ici vaut pour toutes les instances."));
}


// ---------------------------------------------------------------- API
// Le corps d'un calcul ne se fige pas a l'appel de corpsReseau() : il se
// resout DANS api(), apres avoir attendu que l'aplatissement soit a jour.
// Sans cela, un calcul lance juste apres une edition de zone partait sur le
// reseau aplati d'AVANT l'edition -- ou sur rien du tout. Une seule
// resolution, ici, plutot qu'une attente a chacun des dix appels.
const RESEAU = Symbol("reseau");

async function api(chemin, corps) {
  if (corps && corps[RESEAU]) {
    await pret();
    const src = E.reseau ? E.plat : { tables: E.tables, positions: E.positions };
    // les hypotheses de simplification sont celles de l'onglet, et elles seules
    const tables = { ...src.tables, SIMPLIFICATIONS: E.simpAppliquer ? E.simp : {} };
    corps = { tables, positions: src.positions, ...corps };
  }
  const opt = corps
    ? { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(corps) }
    : {};
  const r = await fetch(chemin, opt);
  const txt = await r.text();
  let data;
  try {
    data = JSON.parse(txt);
  } catch {
    throw new Error(txt.slice(0, 300));
  }
  if (!r.ok) throw new Error(typeof data.detail === "string" ? data.detail : JSON.stringify(data.detail));
  return data;
}

const corpsReseau = (extra = {}) => ({ [RESEAU]: true, ...extra });

// --- taches de fond ------------------------------------------------------
// Un calcul long part en TACHE : la requete rend un identifiant, le flux SSE
// donne l'avancement (que le moteur signale lui-meme), et le bouton Arreter
// interrompt au prochain point de controle. La promesse rendue se resout sur
// le resultat, exactement comme api() -- l'appelant ne voit pas la difference.
async function lancerTache(genre, corps) {
  const d = await api("/api/taches/" + genre, corps);
  const zone = $("#tache"), barre = $("#tache-barre"), msg = $("#tache-msg");
  zone.hidden = false;
  barre.value = 0;
  msg.textContent = genre + "...";
  $("#tache-arreter").onclick = () => fetch(`/api/taches/${d.id}`, { method: "DELETE" });
  return new Promise((resoudre, rejeter) => {
    const fin = () => {
      flux.close();
      zone.hidden = true;
    };
    const flux = new EventSource(`/api/taches/${d.id}/flux`);
    flux.onmessage = (ev) => {
      const x = JSON.parse(ev.data);
      barre.value = x.fraction;
      msg.textContent = `${genre} ${Math.round(100 * x.fraction)} %` + (x.message ? ` -- ${x.message}` : "");
      if (x.etat === "fini") {
        fin();
        resoudre(x.resultat);
      } else if (x.etat === "echec") {
        fin();
        rejeter(new Error(x.erreur));
      } else if (x.etat === "annule") {
        fin();
        rejeter(new Error("calcul arrete"));
      }
    };
    flux.onerror = () => {
      // flux coupe (serveur relance ?) : on le dit plutot que d'attendre
      fin();
      rejeter(new Error("suivi du calcul interrompu"));
    };
  });
}

// =====================================================================
// tables : lecture / ecriture
// =====================================================================
const COLS = (table) => (E.defauts ? E.defauts.tables[table] || [] : []);
const lignes = (table) => E.tables[table] || (E.tables[table] = []);

function busIds() {
  return lignes("Y_DER").map((r) => +r[0]);
}

function busIdsDC() {
  return lignes("T_BusDC").map((r) => +r[0]);
}

// Les bus alternatifs et continus PARTAGENT l'espace de numeros : c'est le
// contrat des tables, et le cablage s'appuie dessus pour retrouver un noeud
// sans avoir a savoir son genre. Numeroter les deux familles separement
// ferait deux bus 3, dont l'un ecraserait l'autre a l'assemblage.
function nouveauNumero() {
  const ids = [...busIds(), ...busIdsDC()];
  return ids.length ? Math.max(...ids) + 1 : 1;
}

function nouvelId(table) {
  const ids = lignes(table).map((r) => +r[0]);
  return ids.length ? Math.max(...ids) + 1 : 1;
}

/** Les equipements poses sur un bus, toutes tables confondues. */
function equipementsDuBus(b) {
  const out = [];
  for (const [fam, spec] of Object.entries(E.defauts.equipements)) {
    lignes(spec.table).forEach((r, i) => {
      if (+r[1] !== b) return;
      // un shunt a B < 0 est une REACTANCE : meme table, autre symbole
      const famille = spec.table === "T_Shunts" && +r[2] < 0 ? "Reactance" : fam;
      out.push({ famille, table: spec.table, index: i, couleur: spec.couleur });
    });
  }
  return out;
}

const TABLES_BRANCHE = ["T_Lines_SI", "T_Lines_pu", "T_Transfos"];

// --- disjoncteurs -------------------------------------------------------
// Un disjoncteur ouvert est un couple [table, ID] de E.tables.OUVERTS : les
// tables elles-memes ne bougent pas, et refermer rend le reseau d'origine.
// C'est l'assembleur qui retire ce qui est ouvert et ce qui en devient hors
// tension ; E.disj en garde le compte rendu (/api/disjoncteurs) pour le schema.
const OUVRABLES = [
  "T_Lines_SI", "T_Lines_pu", "T_Transfos", "T_IBs", "T_SMs", "T_GFM1s", "T_GFM2s",
  "T_GFM3s", "T_GFM4s", "T_GFL1s", "T_GFL2s", "T_Loads", "T_Shunts", "Y_DER",
];

function estOuvertId(table, id) {
  return (E.tables.OUVERTS || []).some(([a, b]) => a === table && +b === +id);
}

function estOuvert(table, index) {
  const r = lignes(table)[index];
  return !!r && estOuvertId(table, +r[0]);
}

function horsTension(b) {
  return !!(E.disj && (E.disj.hors_tension || []).includes(+b));
}

function basculerDisjoncteur(table, id, ouvert) {
  const reste = (E.tables.OUVERTS || []).filter(([a, b]) => !(a === table && +b === +id));
  E.tables.OUVERTS = ouvert ? [...reste, [table, +id]] : reste;
  modifie();
}

/** Le symbole d'un disjoncteur ouvert : un carre rouge evide, comme sur un
 *  schema d'exploitation. */
function marqueOuvert(x, y) {
  return el("rect", { x: x - 4.5, y: y - 4.5, width: 9, height: 9, fill: "#fff", stroke: "#bf2617", "stroke-width": 2 });
}

/** Les liens continus a tracer : cables entre bus continus, et le trait qui
 *  relie un convertisseur a son bus continu. Le second n'est pas une branche
 *  au sens du modele -- c'est une LIGNE de T_VSCs, qui declare simplement que
 *  la machine de ce bus alternatif a un port continu -- mais sans le dessiner
 *  le schema montrerait deux moities de reseau sans rapport. */
function branchesDC() {
  const out = [];
  lignes("T_LinesDC").forEach((r, i) =>
    out.push({ table: "T_LinesDC", index: i, de: +r[1], vers: +r[2], genre: "cable" })
  );
  lignes("T_VSCs").forEach((r, i) =>
    out.push({ table: "T_VSCs", index: i, de: +r[1], vers: +r[2], genre: "vsc" })
  );
  return out;
}

/** Les familles qui ont une variante a port continu : ce sont les seules
 *  qu'on peut raccorder. Un GFM1 n'a pas de bus continu a lui, son cote
 *  continu est ideal. */
const FAMILLES_VSC = { T_GFM3s: "GFM3", T_GFM4s: "GFM4", T_GFL2s: "GFL2" };

function convertisseurDuBus(b) {
  for (const [table, fam] of Object.entries(FAMILLES_VSC)) {
    const i = lignes(table).findIndex((r) => +r[1] === b);
    if (i >= 0) return { table, famille: fam, index: i };
  }
  return null;
}

function branches() {
  const out = [];
  for (const t of TABLES_BRANCHE)
    lignes(t).forEach((r, i) =>
      out.push({ table: t, index: i, de: +r[1], vers: +r[2], transfo: t === "T_Transfos" })
    );
  return out;
}

// =====================================================================
// reseau vide et chargement
// =====================================================================
function reseauVide() {
  const t = { f0: 50, Sb: 100e6, nom: "nouveau" };
  for (const nom of Object.keys(E.defauts.tables)) t[nom] = [];
  t.SM_PARAMS = {};
  t.GFM_PARAMS = {};
  t.GFL_PARAMS = {};
  t.PARAMS_PAR_BUS = {};
  return t;
}

/** Les hypotheses de simplification que porte un reseau ouvert deviennent
 *  celles de l'onglet, appliquees : c'est ainsi qu'il a ete enregistre. */
function adopterSimplifications(tables) {
  const s = (tables && tables.SIMPLIFICATIONS) || {};
  E.simp = JSON.parse(JSON.stringify(s));
  E.simpAppliquer = Object.keys(E.simp).length > 0;
  E.simpInvVersion = -1;
  const c = $("#simp-appliquer");
  if (c) c.checked = E.simpAppliquer;
}

function poser(tables, positions, nom) {
  E.multi = null;
  adopterSimplifications(tables);
  delete tables.SIMPLIFICATIONS;
  E.tables = tables;
  E.positions = positions || {};
  if (nom) E.tables.nom = nom;
  $("#nom-reseau").value = E.tables.nom || "";
  E.sel = null;
  E.attente = null;
  // un reseau pose est une entree d'historique a lui seul : ni fondu avec
  // l'edition d'avant, ni avec celle qui suit
  H.t = 0;
  modifie();
  H.t = 0;
}

// ``siVide`` : le chargement du demarrage, qui ne remplit qu'un atelier
// encore vide -- un reseau pose pendant qu'il arrivait (Nouveau, Ouvrir,
// Charger) n'est pas ecrase.
async function chargerCatalogue(nom, siVide = false) {
  etat("chargement de " + nom + "...");
  try {
    const d = await api("/api/reseaux/" + nom);
    if (siVide && E.tables) {
      // le message est peut-etre deja celui de l'action qui a pris la place
      if (I18N.original($("#etat")) === "chargement de " + nom + "...") etat("");
      return;
    }
    if (d.instances && d.instances.length) {
      entrerZones(d);
      etat(`${nom} charge : ${d.tables.Y_DER.length} bus a plat, ${d.instances.length} instances de zone`);
    } else {
      sortirZones();
      poser(d.tables, d.positions);
      etat(nom + " charge : " + lignes("Y_DER").length + " bus");
    }
  } catch (err) {
    etat("echec : " + err.message, true);
  }
}

// =====================================================================
// schema SVG (partage constructeur / load flow)
// =====================================================================
const MARGE = { x: 46, y: 34 };

function versPixels(svg, p) {
  const vb = svg.viewBox.baseVal;
  return [MARGE.x + p[0] * (vb.width - 2 * MARGE.x), MARGE.y + p[1] * (vb.height - 2 * MARGE.y)];
}

function versUnite(svg, px, py) {
  const vb = svg.viewBox.baseVal;
  return [
    Math.min(1, Math.max(0, (px - MARGE.x) / (vb.width - 2 * MARGE.x))),
    Math.min(1, Math.max(0, (py - MARGE.y) / (vb.height - 2 * MARGE.y))),
  ];
}

function pointSVG(svg, ev) {
  // Par la matrice ecran de l'SVG, et non par son rectangle englobant : avec
  // preserveAspectRatio le viewBox est CENTRE dans l'element, avec des bandes
  // de part et d'autre. Supposer qu'il le remplit decale tout ce qu'on pose,
  // d'autant plus que la fenetre s'eloigne du rapport du viewBox.
  const p = svg.createSVGPoint();
  p.x = ev.clientX;
  p.y = ev.clientY;
  const q = p.matrixTransform(svg.getScreenCTM().inverse());
  return [q.x, q.y];
}

function positionBus(b) {
  const p = E.positions[String(b)];
  return p ? [p[0], p[1]] : [0.5, 0.5];
}

const COULEUR_TYPE = { 1: "#1c2330", 2: "#1f4e79", 3: "#8a93a5" };

// =====================================================================
// symboles de schema unifilaire
// =====================================================================
// Un schema unifilaire se lit par ses SYMBOLES. Une machine synchrone est un
// cercle marque d'une sinusoide, un transformateur deux cercles enlaces, un
// convertisseur un carre qui dit ce qu'il est, une charge une fleche qui
// quitte le bus. Des pastilles de couleur, elles, demandaient une legende :
// ces symboles, tout electricien les lit sans.

const ROUGE_SEL = "#bf2617";

// Les symboles "a queue" -- charge, condensateur -- partent du bus dans la
// direction ou on les pose, et tournent avec elle. Les autres sont des
// boites qui restent droites : on ne lit pas "GFM" a l'envers.
const ROTATIFS = new Set(["Charge", "ChargeYDER", "ChargeDC", "Shunt", "ShuntDC", "Reactance"]);

function couleurFamille(fam) {
  const eq = (E.defauts && E.defauts.equipements) || {};
  if (eq[fam]) return eq[fam].couleur;
  if (fam === "ChargeYDER") return (eq.Charge || {}).couleur || "#bf2617";
  if (fam === "ShuntDC" || fam === "Reactance") return (eq.Shunt || {}).couleur || "#2a8fb5";
  if (fam === "ChargeDC") return (eq.Charge || {}).couleur || "#bf2617";
  if (fam === "DCDC") return "#7a4fa3";
  return "#344054";
}

/** Les sources d'energie d'un convertisseur, et la couleur de leur
 *  pictogramme. */
const SOURCES = {
  batterie: { titre: "batterie", couleur: "#2e9e5b" },
  PV: { titre: "photovoltaique", couleur: "#d98a00" },
  DC: { titre: "bus continu (reseau DC)", couleur: "#7a4fa3" },
};

/** La source d'energie d'un GFM, d'un GFL ou d'un DC/DC.
 *
 *  Un convertisseur raccorde a un bus continu (``T_VSCs``) prend la source de
 *  son ile continue : une batterie s'il y en a une, sinon un PV, sinon le bus
 *  continu lui-meme (un reseau DC, ou rien de branche). Un convertisseur a bus
 *  continu IDEAL suit la convention de l'editeur : un GFM est a batterie, un
 *  GFL photovoltaique. */
function sourceEquipement(table, index) {
  const row = lignes(table)[index];
  if (!row) return null;
  if (table === "T_DCDC") return String(row[3]).toLowerCase() === "batterie" ? "batterie" : "PV";
  if (!/^T_GF[ML]\ds$/.test(table)) return null;
  const vsc = lignes("T_VSCs").find((r) => +r[1] === +row[1]);
  if (!vsc) return table.startsWith("T_GFM") ? "batterie" : "PV";
  const ile = ileDC(+vsc[2]);
  const types = lignes("T_DCDC").filter((r) => ile.has(+r[1])).map((r) => String(r[3]).toLowerCase());
  if (types.includes("batterie")) return "batterie";
  if (types.includes("pv")) return "PV";
  return "DC";
}

/** Le pictogramme d'une source, dans une pastille centree en (x, y). */
function pictogrammeSource(source, x, y, encre) {
  const s = SOURCES[source];
  const g = el("g", { class: "pictogramme-source", "data-source": source });
  if (!s) return g;
  const c = s.couleur;
  g.append(el("circle", { cx: x, cy: y, r: 7.5, fill: "#fff", stroke: encre, "stroke-width": 1.3 }));
  if (source === "batterie") {
    // une pile : le corps, la borne, et le plus
    g.append(el("rect", { x: x - 4.8, y: y - 2.8, width: 8.4, height: 5.6, rx: 0.8, fill: c }));
    g.append(el("rect", { x: x + 3.6, y: y - 1.4, width: 1.5, height: 2.8, fill: c }));
    g.append(el("line", { x1: x - 2.6, y1: y, x2: x + 0.6, y2: y, stroke: "#fff", "stroke-width": 1 }));
    g.append(el("line", { x1: x - 1, y1: y - 1.6, x2: x - 1, y2: y + 1.6, stroke: "#fff", "stroke-width": 1 }));
  } else if (source === "PV") {
    // un panneau incline sous un soleil
    g.append(el("circle", { cx: x + 2.6, cy: y - 3.2, r: 1.7, fill: c }));
    g.append(el("polygon", { points: `${x - 5},${y + 4} ${x - 2.5},${y - 1} ${x + 5},${y - 1} ${x + 2.5},${y + 4}`, fill: c }));
    g.append(el("line", { x1: x - 3.75, y1: y + 1.5, x2: x + 3.75, y2: y + 1.5, stroke: "#fff", "stroke-width": 0.7 }));
    g.append(el("line", { x1: x + 0.5, y1: y - 1, x2: x - 1.2, y2: y + 4, stroke: "#fff", "stroke-width": 0.7 }));
  } else {
    // le bus continu : le carre des noeuds continus, et le signe "="
    g.append(el("rect", { x: x - 4, y: y - 4, width: 8, height: 8, fill: c }));
    g.append(el("line", { x1: x - 2.4, y1: y - 1.2, x2: x + 2.4, y2: y - 1.2, stroke: "#fff", "stroke-width": 1 }));
    g.append(el("line", { x1: x - 2.4, y1: y + 1.2, x2: x + 2.4, y2: y + 1.2, stroke: "#fff", "stroke-width": 1 }));
  }
  return g;
}

/** Demi-cote des boites de convertisseur. */
const DEMI_CONV = 15;
const estConvertisseur = (fam) => /^GF[ML]/.test(fam) || fam === "DCDC";

/** La boite d'un convertisseur : un carre coupe par la diagonale qui separe
 *  ses deux faces. ``u`` (unitaire) pointe vers le bus de la face ``avant``
 *  (``~`` pour un onduleur), la face ``arriere`` (``=``) lui est opposee. On
 *  prend la diagonale qui laisse ``u`` dans un seul des deux triangles. */
function boiteConvertisseur(g, encre, ep, u, avant, arriere) {
  const s = DEMI_CONV;
  g.append(el("rect", { x: -s, y: -s, width: 2 * s, height: 2 * s, rx: 2.5, fill: "#fff", stroke: encre, "stroke-width": ep }));
  // "/" separe haut-gauche et bas-droite ; "\" separe haut-droite et bas-gauche
  const slash = Math.abs(u[0] + u[1]) >= Math.abs(u[0] - u[1]);
  if (slash) g.append(el("line", { x1: -s, y1: s, x2: s, y2: -s, stroke: encre, "stroke-width": 1.2 }));
  else g.append(el("line", { x1: -s, y1: -s, x2: s, y2: s, stroke: encre, "stroke-width": 1.2 }));
  // les centres des deux triangles ; la face avant est celui que vise u
  const c1 = slash ? [-s / 2.6, -s / 2.6] : [s / 2.6, -s / 2.6];
  const c2 = [-c1[0], -c1[1]];
  const [cav, car] = c1[0] * u[0] + c1[1] * u[1] > 0 ? [c1, c2] : [c2, c1];
  g.append(el("text", { x: cav[0], y: cav[1] + 4, "text-anchor": "middle", class: "symbole-face", fill: encre }, avant));
  g.append(el("text", { x: car[0], y: car[1] + 4, "text-anchor": "middle", class: "symbole-face", fill: encre }, arriere));
}

/** Le point de la boite d'un convertisseur dans la direction ``d``. */
function borneBoite(c, d) {
  const m = Math.max(Math.abs(d[0]), Math.abs(d[1])) || 1;
  return [c[0] + (d[0] * DEMI_CONV) / m, c[1] + (d[1] * DEMI_CONV) / m];
}

/** Une source d'energie dessinee en entier, centree sur l'origine : une
 *  batterie (deux plaques, le plus et le moins) ou un panneau
 *  photovoltaique sous son soleil. */
function symboleSource(source, encre, sel) {
  const g = el("g", { class: "symbole-source", "data-source": source });
  const c = sel ? ROUGE_SEL : (SOURCES[source] || {}).couleur || encre;
  if (source === "batterie") {
    g.append(el("rect", { x: -13, y: -8, width: 26, height: 16, rx: 2, fill: "#fff", stroke: c, "stroke-width": 1.8 }));
    g.append(el("rect", { x: 13, y: -3.5, width: 3, height: 7, rx: 0.8, fill: c }));
    g.append(el("line", { x1: -4, y1: -5, x2: -4, y2: 5, stroke: c, "stroke-width": 3 }));
    g.append(el("line", { x1: 2, y1: -2.5, x2: 2, y2: 2.5, stroke: c, "stroke-width": 3 }));
    g.append(el("text", { x: -9, y: 3.5, "text-anchor": "middle", class: "symbole-indice", fill: c }, "+"));
    g.append(el("text", { x: 8, y: 3.5, "text-anchor": "middle", class: "symbole-indice", fill: c }, "−"));
  } else if (source === "PV") {
    g.append(el("circle", { cx: 10, cy: -10, r: 3.5, fill: c }));
    for (const a of [0, 1, 2, 3, 4, 5, 6, 7]) {
      const t = (a * Math.PI) / 4;
      g.append(el("line", { x1: 10 + 5 * Math.cos(t), y1: -10 + 5 * Math.sin(t), x2: 10 + 7 * Math.cos(t), y2: -10 + 7 * Math.sin(t), stroke: c, "stroke-width": 1 }));
    }
    g.append(el("polygon", { points: "-14,9 -9,-5 9,-5 4,9", fill: "#fff", stroke: c, "stroke-width": 1.8, "stroke-linejoin": "round" }));
    g.append(el("line", { x1: -11.5, y1: 2, x2: 6.5, y2: 2, stroke: c, "stroke-width": 1 }));
    g.append(el("line", { x1: -3, y1: -5, x2: -6, y2: 9, stroke: c, "stroke-width": 1 }));
    g.append(el("line", { x1: 3, y1: -5, x2: 0, y2: 9, stroke: c, "stroke-width": 1 }));
  }
  return g;
}

/** Le symbole d'une famille, en coordonnees locales centrees sur l'origine.
 *  Les symboles a queue commencent en y = -15 et descendent : pose avec une
 *  rotation, leur queue touche le bus. ``source`` (batterie, PV, DC) ajoute
 *  aux convertisseurs le pictogramme de leur source d'energie. */
function symbole(fam, couleur, sel, source = null, u = null) {
  const encre = sel ? ROUGE_SEL : couleur;
  const ep = sel ? 2.8 : 2;
  const g = el("g");
  const trait = (x1, y1, x2, y2, w = ep) =>
    g.append(el("line", { x1, y1, x2, y2, stroke: encre, "stroke-width": w, "stroke-linecap": "round" }));

  if (fam === "SM") {
    g.append(el("circle", { cx: 0, cy: 0, r: 14, fill: "#fff", stroke: encre, "stroke-width": ep }));
    g.append(el("path", {
      d: "M -8.5 0 C -6.5 -8, -2 -8, 0 0 C 2 8, 6.5 8, 8.5 0",
      fill: "none", stroke: encre, "stroke-width": 1.8, "stroke-linecap": "round",
    }));
  } else if (fam.startsWith("GFM") || fam.startsWith("GFL")) {
    // le convertisseur a deux faces : "~" cote bus alternatif, "=" cote
    // continu ; ``u`` pointe vers le bus alternatif (en haut a gauche dans la
    // palette)
    boiteConvertisseur(g, encre, ep, u || [-Math.SQRT1_2, -Math.SQRT1_2], "~", "=");
    // la famille et sa variante (1 a 4) : un GFM3 et un GFM1 n'ont ni les
    // memes etats ni le meme bus continu, et le schema doit le laisser voir
    if (u) {
      // l'etiquette sur le cote libre (perpendiculaire aux deux faces), du
      // cote du haut ; ancree pour s'eloigner de la boite
      let n = [-u[1], u[0]];
      if (n[1] > 0 || (Math.abs(n[1]) < 1e-9 && n[0] < 0)) n = [-n[0], -n[1]];
      const x = n[0] * 22, y = n[1] * 22 + 3;
      const ancre = n[0] > 0.45 ? "start" : n[0] < -0.45 ? "end" : "middle";
      g.append(el("text", { x, y: ancre === "middle" ? y - 2 : y, "text-anchor": ancre, class: "symbole-indice", fill: encre }, fam));
    }
    if (source) g.append(pictogrammeSource(source, 14, -14, encre));
  } else if (fam === "IB") {
    // le reseau amont : carre hachure, la convention des outils de reseau
    g.append(el("rect", { x: -13, y: -13, width: 26, height: 26, fill: "#fff", stroke: encre, "stroke-width": ep }));
    for (const c of [-18, -9, 0, 9, 18]) {
      const x1 = Math.max(-13, c - 13), x2 = Math.min(13, c + 13);
      trait(x1, x1 - c, x2, x2 - c, 1.2);
    }
  } else if (fam === "Blocs") {
    // un schema-bloc en miniature : deux blocs relies, dans un cadre
    g.append(el("rect", { x: -14, y: -12, width: 28, height: 24, rx: 3, fill: "#fff", stroke: encre, "stroke-width": ep }));
    g.append(el("rect", { x: -10, y: -5, width: 7, height: 10, fill: "none", stroke: encre, "stroke-width": 1.3 }));
    g.append(el("rect", { x: 3, y: -5, width: 7, height: 10, fill: "none", stroke: encre, "stroke-width": 1.3 }));
    trait(-3, 0, 3, 0, 1.3);
    trait(-14, 0, -10, 0, 1.3);
    trait(10, 0, 14, 0, 1.3);
  } else if (fam === "DCDC") {
    boiteConvertisseur(g, encre, ep, u || [-Math.SQRT1_2, -Math.SQRT1_2], "=", "=");
    if (source) g.append(pictogrammeSource(source, 14, -14, encre));
  } else if (fam === "Charge" || fam === "ChargeYDER" || fam === "ChargeDC") {
    trait(0, -15, 0, 3);
    g.append(el("polygon", { points: "-6,1 6,1 0,13", fill: encre, stroke: encre, "stroke-width": 1, "stroke-linejoin": "round" }));
  } else if (fam === "Shunt" || fam === "ShuntDC") {
    trait(0, -15, 0, -4);
    trait(-9, -4, 9, -4, ep + 0.6);
    trait(-9, 1, 9, 1, ep + 0.6);
    trait(0, 1, 0, 6);
    trait(-7, 6, 7, 6, 1.6);
    trait(-4.5, 9.5, 4.5, 9.5, 1.6);
    trait(-2, 13, 2, 13, 1.6);
  } else if (fam === "Reactance") {
    // une bobine : trois spires, puis la terre
    trait(0, -15, 0, -9);
    g.append(el("path", {
      d: "M 0 -9 a 3.2 3.2 0 0 1 0 6.4 a 3.2 3.2 0 0 1 0 6.4 a 3.2 3.2 0 0 1 0 6.4",
      fill: "none", stroke: encre, "stroke-width": ep, "stroke-linecap": "round",
    }));
    trait(0, 10.2, 0, 11);
    trait(-7, 11, 7, 11, 1.6);
    trait(-4.5, 14, 4.5, 14, 1.6);
    trait(-2, 17, 2, 17, 1.6);
  } else {
    g.append(el("circle", { cx: 0, cy: 0, r: 10, fill: "#fff", stroke: encre, "stroke-width": ep }));
  }
  return g;
}

/** Deux cercles enlaces au milieu de la branche. Chacun est d'abord pose
 *  plein, de la couleur du fond, pour masquer la ligne qui le traverse ; les
 *  contours viennent ensuite, sinon le second cercle effacerait le premier. */
function symboleTransfo(a, c, encre, sel) {
  const g = el("g");
  const mx = (a[0] + c[0]) / 2, my = (a[1] + c[1]) / 2;
  const L = Math.hypot(c[0] - a[0], c[1] - a[1]) || 1;
  const ux = (c[0] - a[0]) / L, uy = (c[1] - a[1]) / L;
  for (const k of [-1, 1])
    g.append(el("circle", { cx: mx + k * 5.5 * ux, cy: my + k * 5.5 * uy, r: 8.5, fill: "#fff", stroke: "none" }));
  for (const k of [-1, 1])
    g.append(el("circle", {
      cx: mx + k * 5.5 * ux, cy: my + k * 5.5 * uy, r: 8.5,
      fill: "none", stroke: encre, "stroke-width": sel ? 2.8 : 2,
    }));
  return g;
}

/** Une instance de zone, vue du transport : un bloc qui dit son nom et son
 *  facteur K. Bord pointille : ce n'est pas un equipement, c'est un reseau
 *  entier, qu'on ouvre au double-clic. */
function blocZone(inst, sel) {
  const g = el("g");
  const encre = sel ? ROUGE_SEL : "#7a4fa3";
  g.append(el("rect", { x: -34, y: -15, width: 68, height: 30, rx: 6, fill: "#f6f2fb", stroke: encre, "stroke-width": sel ? 2.6 : 1.8, "stroke-dasharray": "5 3" }));
  const nom = String(inst.zone);
  g.append(el("text", { x: 0, y: -1, "text-anchor": "middle", class: "symbole-texte", fill: encre }, nom.length > 11 ? nom.slice(0, 10) + "\u2026" : nom));
  g.append(el("text", { x: 0, y: 10, "text-anchor": "middle", class: "symbole-indice", fill: encre }, `pref. ${inst.prefixe}${inst.K > 1 ? "  x" + inst.K : ""}`));
  return g;
}

function iconeEquipement(fam) {
  const s = el("svg", { viewBox: "-17 -17 34 34", width: 22, height: 22, class: "icone-symbole" });
  s.append(symbole(fam, couleurFamille(fam), false));
  return s;
}

function iconeTransfo() {
  const s = el("svg", { viewBox: "-17 -17 34 34", width: 22, height: 22, class: "icone-symbole" });
  s.append(el("line", { x1: -17, y1: 0, x2: 17, y2: 0, stroke: "#98a2b3", "stroke-width": 2 }));
  s.append(symboleTransfo([-17, 0], [17, 0], "#344054", false));
  return s;
}

/** Les voisins d'un bus, par les branches alternatives comme continues. */
function voisinsDe(b) {
  const v = [];
  for (const br of [...branches(), ...branchesDC()]) {
    if (br.de === b) v.push(br.vers);
    else if (br.vers === b) v.push(br.de);
  }
  return v;
}

/** La direction ou poser les equipements d'un bus : a l'oppose de ses
 *  voisins, pour que le symbole tombe dans le vide plutot que sur une ligne.
 *  Chaque voisin compte pour un vecteur UNITAIRE -- sinon un voisin lointain
 *  l'emporterait sur un voisin proche, qui est pourtant celui qu'on risque de
 *  recouvrir. Un bus isole pose ses equipements vers le bas. */
function directionLibre(svg, b) {
  const p = versPixels(svg, positionBus(b));
  let dx = 0, dy = 0;
  for (const v of voisinsDe(b)) {
    const q = versPixels(svg, positionBus(v));
    const n = Math.hypot(p[0] - q[0], p[1] - q[1]);
    if (n < 1e-6) continue;
    dx += (p[0] - q[0]) / n;
    dy += (p[1] - q[1]) / n;
  }
  const n = Math.hypot(dx, dy);
  return n < 1e-6 ? [0, 1] : [dx / n, dy / n];
}

function tourner([x, y], a) {
  const c = Math.cos(a), s = Math.sin(a);
  return [c * x - s * y, s * x + c * y];
}

/** Ce qu'il y a a dessiner sur un bus alternatif : les equipements des
 *  tables, et la charge "classique" portee par Y_DER elle-meme -- elle n'est
 *  dans aucune table d'equipement, mais un schema qui ne la montrerait pas
 *  cacherait l'essentiel de la consommation. */
function equipementsDessines(b, r) {
  const eq = equipementsDuBus(b).map((e) => ({ ...e }));
  if (+r[7] || +r[8])
    eq.push({ famille: "ChargeYDER", table: "Y_DER", index: lignes("Y_DER").indexOf(r), couleur: couleurFamille("ChargeYDER") });
  return eq;
}

function texteEquipement(e, b) {
  if (e.famille === "ChargeYDER") {
    const r = lignes("Y_DER").find((x) => +x[0] === b);
    return `charge au bus ${b}` + "\n" + `${fmt(r[7] / 1e6, 2)} MW, ${fmt(r[8] / 1e6, 2)} Mvar`;
  }
  const spec = (E.defauts.equipements || {})[e.famille] || {};
  const row = lignes(e.table)[e.index];
  let t = `${spec.titre || e.famille} #${row[0]} au bus ${b}`;
  if (e.table === "T_Loads") t += "\n" + `${fmt(row[2] / 1e6, 2)} MW, ${fmt(row[3] / 1e6, 2)} Mvar, type ${row[4]}`;
  else if (e.table === "T_Shunts")
    t += "\n" + `B = ${fmt(row[2], 3)} pu (${+row[2] < 0 ? "reactance : branche R-L vers la terre" : "condensateur : capacite du noeud"})`;
  else if (e.table === "T_ShuntsDC") t += "\n" + `C = ${fmt(row[2], 3)} pu`;
  else if (e.table === "T_DCDC") t = `convertisseur DC/DC #${row[0]} (${row[3]}) au bus ${b}` + "\n" + `${fmt(row[4] / 1e6, 2)} MW injectes, ${fmt(row[2] / 1e6, 1)} MW nominaux`;
  else if (e.table === "T_LoadsDC") t = `charge continue #${row[0]} au bus ${b}` + "\n" + `${fmt(row[2] / 1e6, 2)} MW, type ${row[3]}`;
  else t += "\n" + `Sn = ${fmt(row[2] / 1e6, 1)} MVA`;
  const src = sourceEquipement(e.table, e.index);
  if (src && e.table !== "T_DCDC") t += "\n" + `source : ${SOURCES[src].titre}`;
  const vsc = lignes("T_VSCs").find((x) => +x[1] === b);
  if (vsc && e.table !== "T_Loads" && e.table !== "T_Shunts")
    t += "\n" + `raccorde au bus continu ${vsc[2]} (mode ${vsc[3]}` +
      (String(vsc[3]) === "droop" ? `, k = ${fmt(100 * (((E.tables.PARAMS_PAR_BUS || {})[vsc[1]] || {}).kdroop || 0), 3)} %` : "") + ")";
  return t;
}

/** Poser les symboles d'un bus, repartis en eventail autour de la direction
 *  libre. Chaque symbole porte ``data-sel-*`` -- le cliquer selectionne SON
 *  equipement -- et ``data-bus`` -- les outils qui attendent un bus (ligne,
 *  transfo, raccordement) le prennent pour celui de son bus. */
function poserSymboles(g, svg, b, p, eq, interactif, estSel, bornes = null) {
  if (!eq.length) return [0, 1];
  const d = directionLibre(svg, b);
  eq.forEach((e, k) => {
    const dk = tourner(d, (k - (eq.length - 1) / 2) * 0.75);
    const sel = e.famille !== "ChargeYDER" && estSel(e.table, e.index);
    let groupe;
    if (e.famille === "Zone") {
      const cx = p[0] + dk[0] * 54, cy = p[1] + dk[1] * 54;
      g.append(el("line", { x1: p[0], y1: p[1], x2: cx, y2: cy, stroke: sel ? ROUGE_SEL : e.couleur, "stroke-width": 2, "stroke-dasharray": "4 3" }));
      groupe = blocZone(e.inst, sel);
      groupe.setAttribute("transform", `translate(${cx},${cy})`);
      groupe.setAttribute("data-zone", e.inst.zone);
    } else if (estConvertisseur(e.famille)) {
      // la face "~" (ou "=" d'un DC/DC) regarde le bus ; la face continue,
      // a l'oppose, va au bus continu de la liaison ou a la source dessinee
      // raccorde a un bus continu, il se pose sur l'axe des deux bus : la face
      // continue regarde alors le bus continu, et le cable ne traverse pas
      // la boite
      const lien = e.famille !== "DCDC" && FAMILLES_VSC[e.table]
        ? lignes("T_VSCs").find((r) => +r[1] === b) : null;
      let dir = dk;
      if (lien) {
        const q = versPixels(svg, positionBus(+lien[2]));
        const n = Math.hypot(q[0] - p[0], q[1] - p[1]);
        if (n > 1e-6) dir = [(q[0] - p[0]) / n, (q[1] - p[1]) / n];
      }
      const c = [p[0] + dir[0] * 42, p[1] + dir[1] * 42];
      const u = [-dir[0], -dir[1]];
      const ac = borneBoite(c, u), dc = borneBoite(c, dir);
      const encre = sel ? ROUGE_SEL : e.couleur;
      g.append(el("line", { x1: p[0], y1: p[1], x2: ac[0], y2: ac[1], stroke: encre, "stroke-width": 2 }));
      groupe = el("g");
      const boite = symbole(e.famille, e.couleur, sel, null, u);
      boite.setAttribute("transform", `translate(${c[0]},${c[1]})`);
      groupe.append(boite);
      if (lien) {
        if (bornes) bornes[b] = dc;
      } else {
        const src = sourceEquipement(e.table, e.index);
        if (SOURCES[src] && src !== "DC") {
          const s = [dc[0] + dir[0] * 26, dc[1] + dir[1] * 26];
          groupe.append(el("line", { x1: dc[0], y1: dc[1], x2: dc[0] + dir[0] * 12, y2: dc[1] + dir[1] * 12, stroke: "#7a4fa3", "stroke-width": 2 }));
          const icone = symboleSource(src, encre, sel);
          icone.setAttribute("transform", `translate(${s[0]},${s[1]})`);
          groupe.append(icone);
        }
      }
    } else if (ROTATIFS.has(e.famille)) {
      const cx = p[0] + dk[0] * 22, cy = p[1] + dk[1] * 22;
      const angle = (Math.atan2(-dk[0], dk[1]) * 180) / Math.PI;
      groupe = symbole(e.famille, e.couleur, sel, sourceEquipement(e.table, e.index));
      groupe.setAttribute("transform", `translate(${cx},${cy}) rotate(${angle})`);
    } else {
      const cx = p[0] + dk[0] * 38, cy = p[1] + dk[1] * 38;
      g.append(el("line", { x1: p[0], y1: p[1], x2: cx, y2: cy, stroke: sel ? ROUGE_SEL : e.couleur, "stroke-width": 2 }));
      groupe = symbole(e.famille, e.couleur, sel, sourceEquipement(e.table, e.index));
      groupe.setAttribute("transform", `translate(${cx},${cy})`);
    }
    if (e.table !== "__instance" && (estOuvert(e.table, e.index) || horsTension(b))) {
      groupe.setAttribute("opacity", "0.35");
      if (estOuvert(e.table, e.index)) {
        const q = versPixels(svg, positionBus(b));
        g.append(marqueOuvert(q[0] + dk[0] * 12, q[1] + dk[1] * 12));
      }
    }
    groupe.setAttribute("data-sel-table", e.table);
    groupe.setAttribute("data-sel-index", e.index);
    groupe.setAttribute("data-sel-bus", b);
    groupe.setAttribute("data-bus", b);
    if (interactif) groupe.setAttribute("style", "cursor:pointer");
    infobulle(groupe, () =>
      e.famille === "Zone"
        ? `zone ${e.inst.zone}` + "\n" + `prefixe ${e.inst.prefixe}, K = ${e.inst.K}` + "\n" + "double-cliquer pour l'ouvrir"
        : texteEquipement(e, b));
    g.append(groupe);
  });
  return d;
}

/** Dessine le reseau en symboles unifilaires. `couleurBus(b)` colore les
 *  noeuds (le load flow s'en sert pour le module de tension), `etiquette(b)`
 *  ajoute un texte sous chaque bus. */
function dessinerSchema(svg, { couleurBus, etiquette, interactif, transit, couleurBranche } = {}) {
  svg.textContent = "";
  const g = el("g");
  const estSel = (table, index) => !!(E.sel && E.sel.table === table && E.sel.index === index);
  const cliquable = interactif ? "cursor:pointer" : null;

  // Une branche se vise par une bande invisible de 12 px : un trait de 2 px
  // ne se clique pas. La bande porte la selection ; le trait visible, rien.
  const prise = (a, c, table, index, texte) => {
    const n = el("line", {
      x1: a[0], y1: a[1], x2: c[0], y2: c[1], stroke: "transparent", "stroke-width": 12,
      style: cliquable, "data-sel-table": table, "data-sel-index": index,
    });
    infobulle(n, texte);
    return n;
  };

  // --- branches alternatives --------------------------------------------
  for (const b of branches()) {
    const a = versPixels(svg, positionBus(b.de));
    const c = versPixels(svg, positionBus(b.vers));
    const sel = estSel(b.table, b.index);
    const ouverte = estOuvert(b.table, b.index);
    const morte = horsTension(b.de) || horsTension(b.vers);
    const teinte = couleurBranche ? couleurBranche(b) : null;
    const encre = sel ? ROUGE_SEL : ouverte || morte ? "#c9ced6" : teinte ? teinte.couleur : "#98a2b3";
    g.append(el("line", {
      x1: a[0], y1: a[1], x2: c[0], y2: c[1], stroke: encre,
      "stroke-width": sel ? 3.5 : teinte ? teinte.largeur : 2,
      "stroke-dasharray": ouverte ? "6 4" : null,
    }));
    if (b.transfo) g.append(symboleTransfo(a, c, sel ? ROUGE_SEL : "#344054", sel));
    if (ouverte)
      for (const f of [0.18, 0.82])
        g.append(marqueOuvert(a[0] + (c[0] - a[0]) * f, a[1] + (c[1] - a[1]) * f));
    // Le transit : une fleche dans le sens de la puissance active, et sa
    // valeur. Placee aux deux tiers vers le bus qui recoit, pour ne pas
    // tomber sur le symbole du transformateur, qui occupe le milieu.
    const P = transit ? transit(b.de, b.vers) : null;
    if (P !== null && P !== undefined && Math.abs(P) > 1e-6) {
      const [u, w] = P > 0 ? [a, c] : [c, a];
      const L = Math.hypot(w[0] - u[0], w[1] - u[1]) || 1;
      const ux = (w[0] - u[0]) / L, uy = (w[1] - u[1]) / L;
      const t = b.transfo ? 0.78 : 0.68;
      const mx = u[0] + (w[0] - u[0]) * t, my = u[1] + (w[1] - u[1]) * t;
      g.append(el("polygon", {
        points: `${mx + 7 * ux},${my + 7 * uy} ${mx - 5 * ux - 4.5 * uy},${my - 5 * uy + 4.5 * ux} ${mx - 5 * ux + 4.5 * uy},${my - 5 * uy - 4.5 * ux}`,
        fill: "#2a7a4b",
      }));
      g.append(el("text", { x: mx - 10 * uy, y: my + 10 * ux + 3, "text-anchor": "middle", class: "etiquette transit" },
        `${fmt(Math.abs(P) * E.tables.Sb / 1e6, 1)} MW`));
    }
    g.append(prise(a, c, b.table, b.index, () =>
      `${b.transfo ? "transformateur" : "ligne"} #${lignes(b.table)[b.index][0]} : bus ${b.de} - ${b.vers}` +
      (ouverte ? "\n" + "disjoncteur OUVERT" : morte ? "\n" + "hors tension" : "")));
  }

  // --- liens continus ---------------------------------------------------
  // Les cables ici ; les raccordements de convertisseur (T_VSCs) apres les
  // symboles : ils partent de la face continue de la boite, qu'il faut avoir
  // posee pour la connaitre.
  const lienDC = (b, a) => {
    const c = versPixels(svg, positionBus(b.vers));
    const sel = estSel(b.table, b.index);
    g.append(el("line", {
      x1: a[0], y1: a[1], x2: c[0], y2: c[1],
      stroke: sel ? ROUGE_SEL : "#7a4fa3",
      "stroke-width": sel ? 3.5 : b.genre === "cable" ? 2.6 : 2,
    }));
    g.append(prise(a, c, b.table, b.index, () =>
      b.genre === "cable"
        ? `cable continu #${lignes(b.table)[b.index][0]} : bus ${b.de} - ${b.vers}`
        : `VSC #${lignes(b.table)[b.index][0]} : bus alternatif ${b.de} -> bus continu ${b.vers}` +
          "\n" + `mode ${lignes(b.table)[b.index][3]}`));
  };
  for (const b of branchesDC()) {
    if (b.genre === "vsc") continue;
    lienDC(b, versPixels(svg, positionBus(b.de)));
  }
  // --- symboles, SOUS les noeuds : le noeud recouvre le depart du trait ---
  const dirs = {};
  const bornesDC = {};
  const instances = E.reseau && E.vue === "transport" ? E.reseau.instances : [];
  for (const r of lignes("Y_DER")) {
    const b = +r[0];
    const eq = equipementsDessines(b, r);
    // les instances de zone accrochees a ce bus se rangent dans l'eventail,
    // comme un equipement de plus : c'est ce qu'elles sont, vues du transport
    instances.forEach((inst, k) => {
      if (+inst.bus_amont === b)
        eq.push({ famille: "Zone", table: "__instance", index: k, couleur: "#7a4fa3", inst });
    });
    dirs[b] = { d: poserSymboles(g, svg, b, versPixels(svg, positionBus(b)), eq, interactif, estSel, bornesDC), n: eq.length };
  }
  for (const b of branchesDC())
    if (b.genre === "vsc") lienDC(b, bornesDC[b.de] || versPixels(svg, positionBus(b.de)));
  for (const r of lignes("T_BusDC")) {
    const b = +r[0];
    const eq = [];
    lignes("T_ShuntsDC").forEach((sh, i) => {
      if (+sh[1] === b) eq.push({ famille: "ShuntDC", table: "T_ShuntsDC", index: i, couleur: couleurFamille("ShuntDC") });
    });
    lignes("T_DCDC").forEach((r, i) => {
      if (+r[1] === b) eq.push({ famille: "DCDC", table: "T_DCDC", index: i, couleur: couleurFamille("DCDC") });
    });
    lignes("T_LoadsDC").forEach((r, i) => {
      if (+r[1] === b) eq.push({ famille: "ChargeDC", table: "T_LoadsDC", index: i, couleur: couleurFamille("ChargeDC") });
    });
    dirs[b] = { d: poserSymboles(g, svg, b, versPixels(svg, positionBus(b)), eq, interactif, estSel), n: eq.length };
  }

  // Les textes evitent le cote ou partent les symboles : le numero au-dessus
  // du noeud, sauf si les symboles montent ; l'etiquette dessous, ou a droite
  // si ce cote-la est pris.
  const textes = (b, p) => {
    const { d, n } = dirs[b] || { d: [0, 1], n: 0 };
    const monte = n && d[1] < -0.6, descend = n && d[1] > 0.6;
    g.append(el("text", { x: p[0], y: monte ? p[1] + 20 : p[1] - 12, "text-anchor": "middle", class: "etiquette" }, String(b)));
    if (!etiquette) return;
    const attrs = monte || descend
      ? { x: p[0] + 12, y: p[1] + 4, "text-anchor": "start", class: "etiquette" }
      : { x: p[0], y: p[1] + 20, "text-anchor": "middle", class: "etiquette" };
    g.append(el("text", attrs, etiquette(b)));
  };

  // --- bus continus : un carre ------------------------------------------
  for (const r of lignes("T_BusDC")) {
    const b = +r[0];
    const p = versPixels(svg, positionBus(b));
    const sel = E.sel && E.sel.table === "T_BusDC" && +lignes("T_BusDC")[E.sel.index][0] === b;
    const cote = 15;
    const carre = el("rect", {
      x: p[0] - cote / 2, y: p[1] - cote / 2, width: cote, height: cote,
      fill: "#7a4fa3", stroke: sel ? ROUGE_SEL : "#fff", "stroke-width": sel ? 3 : 2,
      class: interactif ? "bus" : null, "data-bus": b, "data-dc": "1",
    });
    g.append(carre);
    infobulle(carre, () =>
      `bus continu ${b} -- ${+r[1] === 1 ? "reference Vdc" : "P impose"}` + "\n" +
      `${fmt(r[2] / 1e3, 1)} kV, ${fmt(r[4] / 1e6, 2)} MW injectes`);
    textes(b, p);
  }

  // --- le PCC d'une zone : la ou le transport la tient --------------------
  const zoneOuverte = E.reseau && E.vue !== "transport" ? E.reseau.zones[E.vue] : null;
  if (zoneOuverte) {
    const q = versPixels(svg, positionBus(+zoneOuverte.bus_pcc));
    g.append(el("circle", { cx: q[0], cy: q[1], r: 13, fill: "none", stroke: "#7a4fa3", "stroke-width": 2, "stroke-dasharray": "3 2" }));
    g.append(el("text", { x: q[0] - 16, y: q[1] + 4, "text-anchor": "end", class: "etiquette", fill: "#7a4fa3" }, "PCC"));
  }

  // --- bus alternatifs : un noeud, colore par son type (ou par le load flow)
  for (const r of lignes("Y_DER")) {
    const b = +r[0];
    const p = versPixels(svg, positionBus(b));
    const sel = E.sel && E.sel.table === "Y_DER" && +lignes("Y_DER")[E.sel.index][0] === b;
    const noeud = el("circle", {
      opacity: horsTension(b) ? 0.35 : null,
      cx: p[0], cy: p[1], r: +r[1] === 1 ? 8 : 7,
      fill: couleurBus ? couleurBus(b) : COULEUR_TYPE[+r[1]],
      stroke: sel ? ROUGE_SEL : "#fff", "stroke-width": sel ? 3 : 2,
      class: interactif ? "bus" : null, "data-bus": b,
    });
    g.append(noeud);
    infobulle(noeud, () => {
      const l = [`bus ${b} -- ${["", "slack", "PV", "PQ"][+r[1]]} -- ${fmt(r[2] / 1e3, 1)} kV`];
      if (etiquette) l.push(etiquette(b));
      return l.join("\n");
    });
    textes(b, p);
  }
  svg.append(g);
}

/** Executer ``fn`` en faisant pointer les tables sur une autre source. */
function avecSource(src, fn) {
  const t = E.tables, p = E.positions, v = E.vue;
  E.tables = src.tables;
  E.positions = src.positions || {};
  E.vue = null; // le reseau aplati n'a ni transport ni zone a dessiner
  try {
    fn();
  } finally {
    E.tables = t;
    E.positions = p;
    E.vue = v;
  }
}

function infobulle(noeud, texte) {
  noeud.addEventListener("mousemove", (ev) => {
    const b = $("#infobulle");
    b.textContent = typeof texte === "function" ? texte() : texte;
    b.style.display = "block";
    b.style.left = ev.clientX + 14 + "px";
    b.style.top = ev.clientY + 14 + "px";
  });
  noeud.addEventListener("mouseleave", () => ($("#infobulle").style.display = "none"));
}

// =====================================================================
// constructeur
// =====================================================================
// La palette : UNE seule liste, alternatif et continu ensemble, rangee par
// role dans le reseau. Un GFM ou un GFL se pose avec sa source d'energie --
// batterie, photovoltaique ou bus continu --, et l'editeur pose avec lui
// toute la chaine continue qu'elle demande (voir poserConvertisseurSource).
const OUTILS_SOURCE = {
  gfm_batterie: { genre: "GFM", source: "batterie" },
  gfm_pv: { genre: "GFM", source: "PV" },
  gfm_dc: { genre: "GFM", source: "DC" },
  gfl_pv: { genre: "GFL", source: "PV" },
  gfl_batterie: { genre: "GFL", source: "batterie" },
  gfl_dc: { genre: "GFL", source: "DC" },
};
const AIDE_SOURCE = {
  batterie: "la batterie tient la tension continue, le convertisseur regle sa puissance",
  PV: "le PV injecte sa puissance (MPPT), le convertisseur tient la tension continue",
  DC: "son bus continu, a relier par un cable continu a un reseau DC ; le convertisseur tient la tension continue tant qu'il est seul",
};
const PALETTE = [
  {
    titre: "Reseau",
    outils: [
      { outil: "bus", titre: "Bus", icone: "Bus" },
      { outil: "bus_dc", titre: "Bus continu", icone: "BusDC" },
      { outil: "ligne", titre: "Ligne", icone: "Ligne" },
      { outil: "transfo", titre: "Transformateur", icone: "Transfo" },
      { outil: "cable_dc", titre: "Cable continu", icone: "Ligne" },
      { outil: "IB", titre: "Reseau amont (bus infini)", icone: "IB" },
    ],
  },
  {
    titre: "Production",
    outils: [
      { outil: "SM", titre: "Machine synchrone", icone: "SM" },
      { outil: "gfm_batterie", titre: "GFM batterie", icone: "GFM", source: "batterie" },
      { outil: "gfm_pv", titre: "GFM photovoltaique", icone: "GFM", source: "PV" },
      { outil: "gfm_dc", titre: "GFM sur bus DC", icone: "GFM", source: "DC" },
      { outil: "gfl_pv", titre: "GFL photovoltaique", icone: "GFL", source: "PV" },
      { outil: "gfl_batterie", titre: "GFL batterie", icone: "GFL", source: "batterie" },
      { outil: "gfl_dc", titre: "GFL sur bus DC", icone: "GFL", source: "DC" },
      { outil: "dcdc", titre: "DC/DC (batterie ou PV)", icone: "DCDC", source: "PV", aide: "sur un bus continu : un PV a puissance constante, que ses proprietes passent en batterie" },
      { outil: "Blocs", titre: "Element a schema-bloc", icone: "Blocs", aide: "Element a schema-bloc : cliquer un bus pour y raccorder le modele courant de l'onglet Blocs (sur son propre bus, derriere un transformateur)" },
    ],
  },
  {
    titre: "Cote continu d'un convertisseur",
    outils: [
      { outil: "cote_batterie", titre: "Batterie", icone: "SourceBatterie", aide: "Batterie : cliquer un GFM/GFL (ou son bus continu) ; elle se branche de son cote continu par un DC/DC, et tient la tension continue par statisme" },
      { outil: "cote_pv", titre: "Photovoltaique", icone: "SourcePV", aide: "Photovoltaique : cliquer un GFM/GFL (ou son bus continu) ; il se branche de son cote continu par un DC/DC et injecte sa puissance (MPPT)" },
    ],
  },
  {
    titre: "Consommation et compensation",
    outils: [
      { outil: "Charge", titre: "Charge", icone: "Charge" },
      { outil: "charge_dc", titre: "Charge continue", icone: "ChargeDC" },
      { outil: "Inductance", titre: "Inductance", icone: "Reactance" },
      { outil: "Condensateur", titre: "Condensateur", icone: "Shunt" },
      { outil: "shunt_dc", titre: "Condensateur continu", icone: "ShuntDC" },
    ],
  },
  {
    titre: "Raccordement",
    outils: [
      { outil: "vsc", titre: "Raccorder un GFM/GFL a un bus continu", icone: "VSC" },
    ],
  },
];
for (const groupe of PALETTE)
  for (const o of groupe.outils)
    if (OUTILS_SOURCE[o.outil])
      o.aide = `${o.titre} : cliquer un bus pour l'y raccorder par un transformateur, ou dans le vide pour le poser seul ; ${AIDE_SOURCE[o.source]}`;
// outil de la palette -> famille de E.defauts.equipements
const FAMILLE_OUTIL = { SM: "SM", GFM: "GFM1", GFL: "GFL1", IB: "IB", Charge: "Charge", Inductance: "Shunt", Condensateur: "Shunt", Blocs: "Blocs" };

function iconePalette(icone, source = null) {
  // la pastille de la source deborde en haut a droite : la vue s'elargit
  const vue = source ? "-17 -24 41 41" : "-17 -17 34 34";
  const s = el("svg", { viewBox: vue, width: source ? 30 : 22, height: source ? 30 : 22, class: "icone-symbole" });
  if (icone === "Bus") s.append(el("circle", { cx: 0, cy: 0, r: 6, fill: "#667085" }));
  else if (icone === "BusDC") s.append(el("rect", { x: -6, y: -6, width: 12, height: 12, fill: "#7a4fa3" }));
  else if (icone === "Ligne") s.append(el("line", { x1: -14, y1: 8, x2: 14, y2: -8, stroke: "#344054", "stroke-width": 2.5, "stroke-linecap": "round" }));
  else if (icone === "Transfo") return iconeTransfo();
  else if (icone === "SourceBatterie") s.append(symboleSource("batterie", "#344054", false));
  else if (icone === "SourcePV") s.append(symboleSource("PV", "#344054", false));
  else if (icone === "VSC") {
    s.append(el("rect", { x: -14, y: -14, width: 28, height: 28, rx: 2.5, fill: "#fff", stroke: "#344054", "stroke-width": 2 }));
    s.append(el("line", { x1: -14, y1: 14, x2: 14, y2: -14, stroke: "#344054", "stroke-width": 1.5 }));
    s.append(el("text", { x: -6, y: -3, "text-anchor": "middle", class: "symbole-indice", fill: "#344054" }, "~"));
    s.append(el("text", { x: 7, y: 10, "text-anchor": "middle", class: "symbole-indice", fill: "#344054" }, "="));
  } else s.append(symbole(icone, couleurFamille(icone), false, source));
  return s;
}

const AIDE = {
  choix_import: "cliquer le bus existant auquel raccorder le reseau importe",
  select: "cliquer pour selectionner, glisser pour deplacer un bus ; maintenir le clic dans le vide pour selectionner un groupe, puis glisser l'un de ses bus",
  bus: "cliquer dans le vide pour poser un bus",
  ligne: "cliquer le bus de depart, puis celui d'arrivee",
  ligne_si: "cliquer le bus de depart, puis celui d'arrivee",
  transfo: "cliquer le primaire, puis le secondaire",
  bus_dc: "cliquer dans le vide pour poser un bus continu",
  cable_dc: "cliquer les deux bus CONTINUS a relier",
  vsc: "cliquer le bus alternatif du convertisseur, puis son bus continu",
  shunt_dc: "cliquer le bus CONTINU qui porte le condensateur",
  dcdc: "cliquer le bus CONTINU du convertisseur DC/DC",
  charge_dc: "cliquer le bus CONTINU de la charge",
  lasso: "tracer un rectangle autour des bus de la zone",
  coller: "cliquer le bus qui recevra l'equipement copie",
};

for (const [nom, o] of Object.entries(OUTILS_SOURCE))
  AIDE[nom] =
    `cliquer un bus pour y raccorder un ${o.genre} (${SOURCES[o.source].titre}) par un transformateur, ` +
    "ou dans le vide pour le poser seul, son bus alternatif a relier";

function outil(nom) {
  E.outil = nom;
  if (nom !== "select") E.multi = null;
  E.attente = null;
  $$(".palette button[data-outil]").forEach((b) => b.classList.toggle("actif", b.dataset.outil === nom));
  $("#canevas").classList.toggle("outil-actif", nom !== "select");
  $("#aide-outil").textContent = AIDE[nom] || "cliquer un bus pour y poser l'equipement";
  dessiner();
}

function dessiner() {
  if (!E.tables) return;
  dessinerCanevas($("#canevas"));
  if (E.attente !== null) {
    const p = versPixels($("#canevas"), positionBus(E.attente));
    $("#canevas").append(el("circle", { cx: p[0], cy: p[1], r: 16, fill: "none", stroke: "#bf2617", "stroke-width": 2 }));
  }
  rendreInspecteur();
}

function selectionner(table, index) {
  E.sel = table === null ? null : { table, index };
  dessiner();
}

// --- gestes -----------------------------------------------------------
function initCanevas() {
  const svg = $("#canevas");

  svg.addEventListener("pointerdown", (ev) => {
    // DOUBLE-CLIC SUR UNE ZONE : il se detecte ici, sur deux pointerdown
    // rapproches, et non par l'evenement dblclick. Le premier clic redessine
    // le schema ; le bloc vise est alors remplace, et le dblclick qui suit
    // part vers un ancetre, jamais vers le bloc -- la zone ne s'ouvrirait pas.
    const zb = ev.target.closest("[data-zone]");
    if (zb && E.reseau) {
      const z = zb.dataset.zone;
      if (E._zoneClic === z && ev.timeStamp - E._tZoneClic < 450) {
        E._zoneClic = null;
        ouvrirVue(z);
        return;
      }
      E._zoneClic = z;
      E._tZoneClic = ev.timeStamp;
    }
    // Un symbole ou une branche se selectionne des le POINTERDOWN, pas au
    // click : le pointerdown redessine le schema, si bien que l'element vise
    // n'existe deja plus quand le click arrive -- et le click ne part alors
    // vers rien. Cliquer un symbole selectionne son equipement, et le glisser
    // deplace son bus.
    if (E.outil === "select") {
      const s = ev.target.closest("[data-sel-table]");
      if (s) {
        // DOUBLE-CLIC : la fiche de l'equipement (equations, schemas). Meme
        // detection que pour les zones, et pour la meme raison.
        const cle = `${s.dataset.selTable}:${s.dataset.selIndex}`;
        if (E._selClic === cle && ev.timeStamp - E._tSelClic < 450) {
          E._selClic = null;
          E.drag = null;
          ouvrirFiche(s.dataset.selTable, +s.dataset.selIndex);
          return;
        }
        E._selClic = cle;
        E._tSelClic = ev.timeStamp;
        const bs = s.dataset.selBus !== undefined ? +s.dataset.selBus : null;
        if (bs !== null && E.multi && E.multi.has(bs)) {
          commencerGroupe(svg, ev);
          return;
        }
        viderGroupe();
        selectionner(s.dataset.selTable, +s.dataset.selIndex);
        if (bs !== null) {
          E.drag = { bus: bs };
          svg.setPointerCapture(ev.pointerId);
        }
        return;
      }
    }
    if (E.outil === "lasso") {
      const [px, py] = pointSVG(svg, ev);
      E.lasso = { x0: px, y0: py, x1: px, y1: py };
      svg.setPointerCapture(ev.pointerId);
      dessinerLasso(svg);
      return;
    }
    const cible = ev.target.closest("[data-bus]");
    const b = cible ? +cible.dataset.bus : null;

    const dc = cible ? cible.dataset.dc === "1" : false;

    if (b !== null && E.outil === "select" && E.multi && E.multi.has(b)) {
      commencerGroupe(svg, ev);
      return;
    }
    if (b !== null && E.outil === "select") {
      viderGroupe();
      const table = dc ? "T_BusDC" : "Y_DER";
      selectionner(table, lignes(table).findIndex((r) => +r[0] === b));
      E.drag = { bus: b };
      svg.setPointerCapture(ev.pointerId);
      return;
    }
    if (b !== null && OUTILS_COTE_DC[E.outil]) {
      // cote continu d'un convertisseur : sur son bus continu, ou sur le
      // convertisseur lui-meme (son bus alternatif), qu'on dote alors d'un
      // bus continu s'il n'en a pas
      poserSourceCoteDC(b, dc, OUTILS_COTE_DC[E.outil]);
      outil("select");
      return;
    }
    if (b !== null && ["shunt_dc", "dcdc", "charge_dc"].includes(E.outil)) {
      if (!dc) {
        etat("cet element se pose sur un bus CONTINU", true);
        return;
      }
      if (E.outil === "shunt_dc") poserShuntDC(b);
      else if (E.outil === "dcdc") poserDCDC(b);
      else poserChargeDC(b);
      outil("select");
      return;
    }
    if (b !== null && E.outil === "cable_dc") {
      if (!dc) {
        etat("un cable continu relie deux bus CONTINUS", true);
        return;
      }
      if (E.attente === null) {
        E.attente = b;
        dessiner();
      } else if (E.attente !== b) {
        poserCableDC(E.attente, b);
        E.attente = null;
        outil("select");
      }
      return;
    }
    if (b !== null && E.outil === "vsc") {
      // On attend un bus ALTERNATIF puis un bus CONTINU, dans cet ordre :
      // c'est le sens de lecture de T_VSCs (BusAC, BusDC).
      if (E.attente === null) {
        if (dc) {
          etat("cliquer d'abord le bus alternatif du convertisseur", true);
          return;
        }
        E.attente = b;
        dessiner();
      } else {
        if (!dc) {
          etat("cliquer maintenant le bus CONTINU de raccordement", true);
          return;
        }
        poserVSC(E.attente, b);
        E.attente = null;
        outil("select");
      }
      return;
    }
    if (b !== null && dc) return;  // les outils alternatifs ignorent le continu
    if (b !== null && ["ligne", "ligne_si", "transfo"].includes(E.outil)) {
      if (E.attente === null) {
        E.attente = b;
        dessiner();
      } else if (E.attente !== b) {
        poserBranche(E.attente, b);
        E.attente = null;
        outil("select");
      }
      return;
    }
    if (b !== null && OUTILS_SOURCE[E.outil]) {
      const o = OUTILS_SOURCE[E.outil];
      poserConvertisseurSource(o.genre, o.source, b);
      outil("select");
      return;
    }
    if (b !== null && FAMILLE_OUTIL[E.outil]) {
      poserEquipement(FAMILLE_OUTIL[E.outil], b, E.outil);
      outil("select");
      return;
    }
    if (b !== null && E.outil === "choix_import") {
      $("#imp-bus-existant").value = String(b);
      outil("select");
      return;
    }
    if (b !== null && E.outil === "coller") {
      collerSur(b);
      outil("select");
      return;
    }
    if (b === null && E.outil === "bus") {
      const [px, py] = pointSVG(svg, ev);
      poserBus(versUnite(svg, px, py));
      outil("select");
      return;
    }
    if (b === null && E.outil === "bus_dc") {
      const [px, py] = pointSVG(svg, ev);
      poserBusDC(versUnite(svg, px, py));
      outil("select");
      return;
    }
    if (b === null && OUTILS_SOURCE[E.outil]) {
      const [px, py] = pointSVG(svg, ev);
      const o = OUTILS_SOURCE[E.outil];
      poserConvertisseurSource(o.genre, o.source, null, versUnite(svg, px, py));
      outil("select");
      return;
    }
    if (b === null && E.outil === "select" && !ev.target.closest("[data-zone]")) {
      // maintenir le clic dans le vide : un rectangle de selection
      const [px, py] = pointSVG(svg, ev);
      E.cadreSel = { x0: px, y0: py, x1: px, y1: py };
      svg.setPointerCapture(ev.pointerId);
      return;
    }
    if (b === null) selectionner(null);
  });

  svg.addEventListener("pointermove", (ev) => {
    if (E.cadreSel) {
      const [px, py] = pointSVG(svg, ev);
      E.cadreSel.x1 = px;
      E.cadreSel.y1 = py;
      dessinerCanevas(svg);
      const l = E.cadreSel;
      svg.append(el("rect", {
        x: Math.min(l.x0, l.x1), y: Math.min(l.y0, l.y1),
        width: Math.abs(l.x1 - l.x0), height: Math.abs(l.y1 - l.y0), class: "cadre-sel",
      }));
      for (const b of busDansRectangle(svg, l, true)) halo(svg, b);
      return;
    }
    if (E.drag && E.drag.groupe) {
      const [px, py] = pointSVG(svg, ev);
      const [ux, uy] = versUnite(svg, px, py);
      const g = E.drag;
      // le deplacement est borne pour que le groupe entier reste dans le cadre
      const dx = Math.min(Math.max(ux - g.depart[0], g.dxMin), g.dxMax);
      const dy = Math.min(Math.max(uy - g.depart[1], g.dyMin), g.dyMax);
      for (const [b, p0] of Object.entries(g.origines)) E.positions[b] = [p0[0] + dx, p0[1] + dy];
      dessinerCanevas(svg);
      return;
    }
    if (E.lasso) {
      const [px, py] = pointSVG(svg, ev);
      E.lasso.x1 = px;
      E.lasso.y1 = py;
      dessinerLasso(svg);
      return;
    }
    if (!E.drag) return;
    const [px, py] = pointSVG(svg, ev);
    E.positions[String(E.drag.bus)] = versUnite(svg, px, py);
    dessinerCanevas(svg);
  });

  const fin = () => {
    if (E.cadreSel) {
      const l = E.cadreSel;
      E.cadreSel = null;
      const petit = Math.abs(l.x1 - l.x0) < 6 && Math.abs(l.y1 - l.y0) < 6;
      const bus = petit ? [] : busDansRectangle(svg, l, true);
      E.multi = bus.length ? new Set(bus) : null;
      selectionner(null);
      if (bus.length)
        etat(`${bus.length} bus selectionnes (${bus.join(", ")}) : glisser l'un d'eux deplace le groupe ; Echap ou clic dans le vide pour deselectionner`);
      return;
    }
    if (E.lasso) {
      const l = E.lasso;
      E.lasso = null;
      finirLasso(svg, l);
      return;
    }
    if (E.drag) {
      E.drag = null;
      E.version += 1; // la position n'affecte aucun calcul : pas d'invalidation
      sauverBrouillon();
      memoriser();
    }
  };
  svg.addEventListener("pointerup", fin);
  svg.addEventListener("pointercancel", fin);
}

// --- selection multiple : un groupe de bus qu'on deplace ensemble -----------
function viderGroupe() {
  E.multi = null;
}

function halo(svg, b) {
  const p = versPixels(svg, positionBus(b));
  svg.append(el("circle", { cx: p[0], cy: p[1], r: 12, class: "bus-groupe" }));
}

/** Le schema, et le halo des bus du groupe s'il y en a un. */
function dessinerCanevas(svg) {
  dessinerSchema(svg, { interactif: true });
  if (E.multi) for (const b of E.multi) halo(svg, b);
}

function commencerGroupe(svg, ev) {
  const [px, py] = pointSVG(svg, ev);
  const origines = {};
  for (const b of E.multi) origines[String(b)] = [...positionBus(b)];
  const xs = Object.values(origines).map((p) => p[0]), ys = Object.values(origines).map((p) => p[1]);
  E.drag = {
    groupe: true,
    depart: versUnite(svg, px, py),
    origines,
    dxMin: -Math.min(...xs), dxMax: 1 - Math.max(...xs),
    dyMin: -Math.min(...ys), dyMax: 1 - Math.max(...ys),
  };
  svg.setPointerCapture(ev.pointerId);
}

// --- definir une zone : un rectangle autour des bus ------------------------
function dessinerLasso(svg) {
  dessinerSchema(svg, { interactif: true });
  const l = E.lasso;
  if (!l) return;
  const x = Math.min(l.x0, l.x1), y = Math.min(l.y0, l.y1);
  svg.append(el("rect", { x, y, width: Math.abs(l.x1 - l.x0), height: Math.abs(l.y1 - l.y0), class: "lasso" }));
  for (const b of busDansRectangle(svg, l)) {
    const p = versPixels(svg, positionBus(b));
    svg.append(el("circle", { cx: p[0], cy: p[1], r: 11, class: "bus-lasso" }));
  }
}

function busDansRectangle(svg, l, avecContinu = false) {
  const x0 = Math.min(l.x0, l.x1), x1 = Math.max(l.x0, l.x1);
  const y0 = Math.min(l.y0, l.y1), y1 = Math.max(l.y0, l.y1);
  return (avecContinu ? [...busIds(), ...busIdsDC()] : busIds()).filter((b) => {
    const p = versPixels(svg, positionBus(b));
    return p[0] >= x0 && p[0] <= x1 && p[1] >= y0 && p[1] <= y1;
  });
}

/** Les bus entoures deviennent une ZONE : le serveur la decoupe (une seule
 *  branche vers le reste, pas de slack), l'instancie a la meme place, et le
 *  constructeur passe en mode zones -- la vue Transport / distribution s'en
 *  sert directement. */
async function finirLasso(svg, l) {
  const bus = busDansRectangle(svg, l);
  outil("select");
  if (!bus.length) return;
  if (E.reseau && E.vue !== "transport") {
    etat("une zone se definit depuis le transport, pas depuis une autre zone", true);
    return;
  }
  const nom = prompt(`Nom de la zone (${bus.length} bus : ${bus.join(", ")})`, `zone${Object.keys((E.reseau || {}).zones || {}).length + 1}`);
  if (!nom) return;
  try {
    const reseau = E.reseau ? corpsStructure() : { tables: E.tables, positions: E.positions };
    const d = await api("/api/zones/definir", { reseau, bus, nom });
    entrerZones(d);
    etat(`zone ${nom} creee : ${bus.length} bus, raccordee au bus ${d.instances[d.instances.length - 1].bus_amont}`);
  } catch (err) {
    etat(err.message, true);
  }
}

// --- importer un reseau -------------------------------------------------------
// Le reseau importe arrive renumerote a droite du schema ; sa reference passe
// au reseau existant, et on le raccorde (transfo, ligne, fusion) en un geste.
function ouvrirImport() {
  if (E.reseau && E.vue !== "transport") ouvrirVue("transport");
  const src = $("#imp-source");
  src.textContent = "";
  src.append(el("option", { value: "@fichier" }, "fichier .json..."));
  for (const o of $("#catalogue").options) if (o.value) src.append(el("option", { value: o.value }, I18N.original(o)));
  if (src.options.length > 1) src.selectedIndex = 1;
  $("#panneau-import").hidden = false;
  $("#imp-messages").textContent = "";
  remplirBusExistants();
  changerSourceImport();
  majModeImport();
}

function remplirBusExistants() {
  const s = $("#imp-bus-existant");
  const garde = s.value;
  s.textContent = "";
  for (const r of lignes("Y_DER")) s.append(el("option", { value: r[0] }, `${r[0]} (${fmt(r[2] / 1e3, 1)} kV)`));
  if (garde) s.value = garde;
}

async function changerSourceImport() {
  const v = $("#imp-source").value;
  $("#imp-c-fichier").hidden = v !== "@fichier";
  E.importSource = null;
  $("#imp-bus-importe").textContent = "";
  try {
    let d;
    if (v === "@fichier") {
      const f = $("#imp-fichier").files[0];
      if (!f) return;
      const brut = JSON.parse(await f.text());
      d = await api("/api/fichiers/lire", { reseau: brut });
      E.importSource = brut;
    } else {
      d = await api("/api/reseaux/" + v);
      E.importSource = { catalogue: v };
    }
    // par defaut : le bus du reseau amont (bus infini), sinon le slack
    const ib = new Set((d.tables.T_IBs || []).map((r) => +r[1]));
    const s = $("#imp-bus-importe");
    for (const r of d.tables.Y_DER) {
      const role = ib.has(+r[0]) ? " -- bus infini" : +r[1] === 1 ? " -- slack" : "";
      s.append(el("option", { value: r[0] }, `${r[0]} (${fmt(r[2] / 1e3, 1)} kV)${role}`));
    }
    const defaut = d.tables.Y_DER.find((r) => ib.has(+r[0])) || d.tables.Y_DER.find((r) => +r[1] === 1);
    if (defaut) s.value = defaut[0];
  } catch (err) {
    $("#imp-messages").textContent = "";
    $("#imp-messages").append(el("div", { class: "pb erreur" }, err.message));
  }
}

function majModeImport() {
  const aucun = $("#imp-mode").value === "aucun";
  $$(".imp-bus").forEach((n) => { n.hidden = aucun; });
}

async function importerReseau() {
  if (!E.importSource) {
    etat("choisir d'abord le reseau a importer", true);
    return;
  }
  const mode = $("#imp-mode").value;
  const raccord = mode === "aucun" ? null : {
    mode,
    bus_importe: +$("#imp-bus-importe").value,
    bus_existant: +$("#imp-bus-existant").value,
  };
  const reseau = E.reseau ? corpsStructure() : { nom: $("#nom-reseau").value, tables: E.tables, positions: E.positions };
  const n = $("#imp-messages");
  n.textContent = "";
  try {
    const d = await api("/api/reseau/importer", {
      reseau,
      source: E.importSource,
      raccord,
      garder_slack: $("#imp-garder-slack").checked,
      garder_ib: $("#imp-garder-ib").checked,
    });
    if (d.instances && d.instances.length) entrerZones(d);
    else {
      sortirZones();
      poser(d.tables, d.positions, d.nom);
    }
    for (const m of d.import.messages) n.append(el("div", { class: "pb ok" }, m));
    remplirBusExistants();
    etat(`reseau importe : bus renumerotes a partir de ${d.import.prefixe + 1}`);
  } catch (err) {
    n.append(el("div", { class: "pb erreur" }, err.message));
  }
}

// --- copier / coller ---------------------------------------------------------
// On copie un EQUIPEMENT avec ses parametres (la ligne de sa table, ses
// reglages par bus, et le transfo qui le raccorde) ; le coller sur un bus en
// pose un double, sur son propre bus terminal comme a la creation.
const TABLES_SUR_BUS = ["T_Loads", "T_Shunts", "T_DCDC", "T_LoadsDC", "T_ShuntsDC"];

function copier() {
  if (!E.sel) return;
  const { table, index } = E.sel;
  const ligne = lignes(table)[index];
  const fam = Object.entries(E.defauts.equipements).find(([, sp]) => sp.table === table);
  if (!ligne || (!fam && !TABLES_SUR_BUS.includes(table))) {
    etat("seuls les equipements se copient (machines, convertisseurs, charges, shunts, DC/DC)", true);
    return;
  }
  const bus = +ligne[1];
  const raccord = lignes("T_Transfos").find((r) => +r[2] === bus || +r[1] === bus);
  E.presse = {
    table,
    famille: fam ? fam[0] : null,
    ligne: JSON.parse(JSON.stringify(ligne)),
    reglages: JSON.parse(JSON.stringify((E.tables.PARAMS_PAR_BUS || {})[String(bus)] || null)),
    transfo: raccord && !TABLES_SUR_BUS.includes(table) ? JSON.parse(JSON.stringify(raccord)) : null,
    // la consigne de l'unite (V, Pg, Qg) est sur la ligne Y_DER de son bus
    consigne: JSON.parse(JSON.stringify(lignes("Y_DER").find((r) => +r[0] === bus) || null)),
  };
  etat(`${table} #${ligne[0]} copie : "Coller" puis cliquer un bus`);
}

function coller() {
  if (!E.presse) {
    etat("rien a coller : selectionner un equipement puis Copier", true);
    return;
  }
  outil("coller");
}

function collerSur(b) {
  const P = E.presse;
  if (!P) return;
  const dc = lignes("T_BusDC").some((r) => +r[0] === b);
  const veutDC = ["T_DCDC", "T_LoadsDC", "T_ShuntsDC"].includes(P.table);
  if (dc !== veutDC) {
    etat(`cet equipement se colle sur un bus ${veutDC ? "continu" : "alternatif"}`, true);
    return;
  }
  if (TABLES_SUR_BUS.includes(P.table)) {
    const l = JSON.parse(JSON.stringify(P.ligne));
    l[0] = nouvelId(P.table);
    l[1] = b;
    lignes(P.table).push(l);
    selectionner(P.table, lignes(P.table).length - 1);
    modifie();
    return;
  }
  // une unite : son bus terminal et son transfo, comme a la creation...
  poserEquipement(P.famille, b);
  const t = lignes(P.table);
  const neuf = t[t.length - 1];
  const terminal = +neuf[1];
  // ...puis les valeurs de l'original, colonne par colonne (sauf ID et bus)
  for (let j = 2; j < P.ligne.length; j++) neuf[j] = P.ligne[j];
  if (P.reglages) (E.tables.PARAMS_PAR_BUS = E.tables.PARAMS_PAR_BUS || {})[String(terminal)] = JSON.parse(JSON.stringify(P.reglages));
  if (P.transfo) {
    const tr = lignes("T_Transfos");
    const r = tr[tr.length - 1];
    r[3] = P.transfo[3]; // Sn
    r[6] = P.transfo[6]; // R
    r[7] = P.transfo[7]; // X
    // tension du bus terminal : celle du cote unite de l'original
    const cote = +P.transfo[2] === +P.ligne[1] ? P.transfo[5] : P.transfo[4];
    r[5] = cote;
    const y = lignes("Y_DER").find((x) => +x[0] === terminal);
    if (y) y[2] = cote;
  }
  const y = lignes("Y_DER").find((x) => +x[0] === terminal);
  if (y && P.consigne && !TABLES_SUR_BUS.includes(P.table)) {
    y[3] = P.consigne[3]; // V
    y[5] = P.consigne[5]; // Pg
    y[6] = P.consigne[6]; // Qg
  }
  modifie();
}

function poserBus(pos) {
  // Un bus nu est toujours PQ : c'est l'EQUIPEMENT qui designe le slack, quand
  // on en pose le premier. Faire du premier bus dessine le slack donnerait un
  // reseau dont la reference d'angle n'est portee par aucune source, et le
  // modele d'etat n'aurait alors pas de repere.
  const b = nouveauNumero();
  const un = lignes("Y_DER").length ? lignes("Y_DER")[0][2] : 20e3;
  lignes("Y_DER").push([b, 3, un, 1.0, 0, 0, 0, 0, 0]);
  E.positions[String(b)] = pos;
  selectionner("Y_DER", lignes("Y_DER").length - 1);
  modifie();
}

/** Poser un equipement : il cree SON bus terminal et son transfo de
 *  raccordement, comme le veut la convention des reseaux de test (une unite
 *  vit sur son propre bus, derriere un transfo). Les charges et les gradins,
 *  eux, se posent directement sur le bus vise. */
function poserEquipement(famille, bus, outilNom = famille) {
  const spec = E.defauts.equipements[famille];
  const table = spec.table;

  if (table === "T_Shunts") {
    // meme table pour les deux : B > 0 condensateur, B < 0 inductance
    lignes(table).push([nouvelId(table), bus, outilNom === "Inductance" ? -0.1 : 0.1]);
  } else if (table === "T_Loads") {
    lignes(table).push([nouvelId(table), bus, 1e6, 0.2e6, "P"]);
  } else {
    const hote = lignes("Y_DER").find((r) => +r[0] === bus);
    const b = nouveauNumero();
    const un = famille === "IB" ? hote[2] : 20e3;
    // La reference d'angle se compte PAR ILE, pas pour tout le reseau : une
    // liaison a courant continu transmet de la puissance mais aucun angle,
    // si bien que chaque cote garde la sienne. Sur un HVDC a deux terminaux,
    // la regle "le premier equipement pose devient LA reference" laissait le
    // second cote sans repere, et le load flow le refusait trois gestes plus
    // tard sans qu'on sache pourquoi.
    // Une ZONE n'a pas de reference d'angle : c'est le transport qui la lui
    // donne, par son bus amont. Y poser un slack la ferait refuser a
    // l'aplatissement -- on ne le fait donc jamais dans une vue de zone.
    const enZone = E.reseau && E.vue !== "transport";
    // Par convention, un GFL est un PV (photovoltaique) : il suit, a sa
    // puissance MPPT, sans tenir la tension -- un bus PQ, jamais le slack.
    // Un GFM est a batterie et une SM tient sa tension : bus PV, ou slack
    // s'il est le premier formeur de son ile.
    // un element a schema-bloc qui ecrit son angle et sa vitesse (sorties
    // theta et w) est un FORMEUR, comme un GFM ; sinon il suit, comme un GFL
    const suiveur = famille.startsWith("GFL") || (famille === "Blocs" && !modeleBlocsFormeur(modeleBlocsCourant()));
    const slack = suiveur ? 3 : enZone ? 2 : ileDuBus(bus).some((x) => {
      const r = lignes("Y_DER").find((y) => +y[0] === x);
      return r && +r[1] === 1;
    })
      ? 2
      : 1;
    lignes("Y_DER").push([b, slack, un, 1.0, 0, 0, 0, 0, 0]);
    const p = positionBus(bus);
    E.positions[String(b)] = placeLibre(p, b);
    if (famille === "IB") lignes(table).push([nouvelId(table), b, E.tables.Sb, 1.0, 0, 0.001, 0.05]);
    else if (famille === "SM") lignes(table).push([nouvelId(table), b, E.tables.Sb]);
    else if (famille === "Blocs") lignes(table).push([nouvelId(table), b, modeleBlocsCourant(), E.tables.Sb]);
    else lignes(table).push([nouvelId(table), b, E.tables.Sb, 40e3]);
    lignes("T_Transfos").push([nouvelId("T_Transfos"), bus, b, E.tables.Sb, hote[2], un, 0.0, 0.05]);
  }
  selectionner(table, lignes(table).length - 1);
  modifie();
}

/** Les bus alternatifs relies a ``b`` par des branches, ``b`` compris. */
function ileDuBus(b) {
  const voisins = {};
  for (const x of busIds()) voisins[x] = [];
  for (const br of branches()) {
    if (voisins[br.de] && voisins[br.vers]) {
      voisins[br.de].push(br.vers);
      voisins[br.vers].push(br.de);
    }
  }
  const vus = new Set();
  const pile = [b];
  while (pile.length) {
    const x = pile.pop();
    if (vus.has(x)) continue;
    vus.add(x);
    for (const v of voisins[x] || []) if (!vus.has(v)) pile.push(v);
  }
  return [...vus];
}

function poserBusDC(pos) {
  // Le premier bus continu pose est la REFERENCE de tension continue de son
  // ile ; les suivants sont a puissance imposee. C'est le decalque de la
  // regle alternative, et sans reference le load flow continu n'a pas de
  // point de depart.
  const b = nouveauNumero();
  const type = lignes("T_BusDC").length ? 2 : 1;
  const un = lignes("T_BusDC").length ? lignes("T_BusDC")[0][2] : 40e3;
  lignes("T_BusDC").push([b, type, un, 1.0, 0]);
  E.positions[String(b)] = pos;
  selectionner("T_BusDC", lignes("T_BusDC").length - 1);
  modifie();
}

function poserCableDC(de, vers) {
  lignes("T_LinesDC").push([nouvelId("T_LinesDC"), de, vers, 100.0, 0.05, 0.2, 0.2]);
  const msg = accorderIleDC(de);
  selectionner("T_LinesDC", lignes("T_LinesDC").length - 1);
  modifie();
  if (msg) etat(msg);
}

/** Un GFM ou un GFL AVEC SA SOURCE D'ENERGIE : batterie, photovoltaique ou
 *  bus continu. L'editeur pose toute la chaine que la source demande :
 *
 *    bus alternatif -- convertisseur -- bus continu (+ condensateur)
 *                                         \-- DC/DC batterie ou PV
 *
 *  Cliquer un bus (``hote``) y raccorde le bus alternatif du convertisseur par
 *  un transformateur ; cliquer dans le vide (``pos``) le laisse ouvert, a
 *  relier par une ligne ou un transfo. Avec la source « bus DC », rien n'est
 *  pose cote continu : on y tire un cable continu vers un reseau DC.
 *
 *  La commande suit la source (voir :func:`accorderIleDC`) : batterie ->
 *  la batterie tient la tension continue et le convertisseur sa puissance
 *  (GFM : PV ; GFL : PQ) ; PV ou bus DC -> le convertisseur tient la tension
 *  continue (VdcV / VdcQ).
 *
 *  Le bus continu recoit d'office son condensateur : sans lui, la boucle de
 *  courant du DC/DC et la regulation du GFM se renvoient la tension du bus,
 *  et le mode bus continu / courant DC/DC part instable (Re = +57 1/s pour un
 *  GFM a batterie qui exporte 10 MW ; stable des 0,2 pu).
 *
 *  UNE EXCEPTION, la reference d'angle. Un convertisseur raccorde au continu
 *  tient deux grandeurs, et sur le slack les deux puissances sont libres : il
 *  ne peut pas l'etre. Un GFM batterie qui serait le premier formeur de son
 *  ile se pose donc a bus continu IDEAL (GFM1, la batterie par convention de
 *  l'editeur), qui, lui, peut tenir l'angle. */
function poserConvertisseurSource(genre, source, hote, pos = null) {
  const gfm = genre === "GFM";
  const Sb = E.tables.Sb;
  const enZone = E.reseau && E.vue !== "transport";
  const sansReference = hote !== null && !enZone && !ileDuBus(hote).some((x) => {
    const r = lignes("Y_DER").find((y) => +y[0] === x);
    return r && +r[1] === 1;
  });
  if (gfm && source === "batterie" && sansReference) {
    poserEquipement("GFM1", hote, "GFM");
    etat(
      `GFM batterie pose au bus ${lignes("Y_DER").at(-1)[0]} : premier formeur de son ile, il en tient ` +
      "la reference d'angle (slack), d'ou un bus continu ideal -- la batterie"
    );
    return;
  }
  const bAC = nouveauNumero();
  // un GFM tient sa tension (bus PV), un GFL non (bus PQ) ; jamais le slack.
  // Derriere une batterie, le convertisseur regle sa puissance : 10 % de Sb.
  const pg = source === "batterie" ? Sb * 0.1 : 0;
  lignes("Y_DER").push([bAC, gfm ? 2 : 3, 20e3, 1.0, 0, pg, 0, 0, 0]);
  if (hote !== null) {
    const h = lignes("Y_DER").find((r) => +r[0] === hote);
    E.positions[String(bAC)] = placeLibre(positionBus(hote), bAC);
    lignes("T_Transfos").push([nouvelId("T_Transfos"), hote, bAC, Sb, h[2], 20e3, 0.0, 0.05]);
  } else E.positions[String(bAC)] = pos;
  const bDC = nouveauNumero();
  const unDC = lignes("T_BusDC").length ? lignes("T_BusDC")[0][2] : 40e3;
  lignes("T_BusDC").push([bDC, 1, unDC, 1.0, 0]);
  E.positions[String(bDC)] = placeLibre(E.positions[String(bAC)], bDC);
  lignes("T_ShuntsDC").push([nouvelId("T_ShuntsDC"), bDC, 0.2]);
  const table = gfm ? "T_GFM3s" : "T_GFL2s";
  lignes(table).push([nouvelId(table), bAC, Sb, unDC]);
  lignes("T_VSCs").push([nouvelId("T_VSCs"), bAC, bDC, gfm ? "VdcV" : "VdcQ"]);
  // la source : un DC/DC de Sb/2 nominaux, qui injecte 10 % de Sb (le PV a
  // sa puissance MPPT ; la batterie, a statisme, se cale sur la tension)
  if (source !== "DC") lignes("T_DCDC").push([nouvelId("T_DCDC"), bDC, Sb * 0.5, source === "PV" ? "PV" : "batterie", Sb * 0.1]);
  accorderIleDC(bDC);
  selectionner(table, lignes(table).length - 1);
  modifie();
  const ac = hote !== null ? `raccorde au bus ${hote} par un transformateur` : "a relier au reseau par une ligne ou un transfo";
  const dc = source === "DC" ? `bus continu ${bDC}, a relier a un reseau DC par un cable continu` : `bus continu ${bDC} avec son DC/DC ${SOURCES[source].titre}`;
  const mode = lignes("T_VSCs").at(-1)[3];
  let msg = `${genre} ${SOURCES[source].titre} pose : bus alternatif ${bAC} ${ac} ; ${dc} ; commande ${mode}`;
  if (sansReference) msg += ". Attention : son ile n'a pas de reference d'angle -- y poser un reseau amont, une SM ou un GFM batterie";
  etat(msg, sansReference);
}

/** Les bus continus relies a ``b`` par des cables, ``b`` compris. */
function ileDC(b) {
  const vus = new Set();
  const pile = [b];
  while (pile.length) {
    const x = pile.pop();
    if (vus.has(x)) continue;
    vus.add(x);
    for (const r of lignes("T_LinesDC")) {
      if (+r[1] === x && !vus.has(+r[2])) pile.push(+r[2]);
      if (+r[2] === x && !vus.has(+r[1])) pile.push(+r[1]);
    }
  }
  return vus;
}

/** Accorder la commande des convertisseurs d'une ile continue a ce qui y est
 *  branche. Une ile continue veut UNE reference de tension :
 *
 *  * un seul convertisseur et une batterie : la batterie tient la tension
 *    continue, le convertisseur regle sa puissance (GFM : PV, statisme P-f ;
 *    GFL : PQ) -- c'est le stockage derriere un onduleur ;
 *  * sinon (PV seul, rien, ou un cable vers d'autres convertisseurs) : un
 *    convertisseur, et un seul, tient la tension continue (VdcV / VdcQ) ; un
 *    PV a puissance constante ne sait pas le faire.
 *
 *  Le bus continu de la reference devient le type 1. Renvoie un message
 *  quand une commande a change, ``null`` sinon. */
function accorderIleDC(busDC) {
  const ile = ileDC(busDC);
  const vscs = lignes("T_VSCs").filter((r) => ile.has(+r[2]));
  if (!vscs.length) return null;
  const gfl = (r) => {
    const c = convertisseurDuBus(+r[1]);
    return c ? c.famille.startsWith("GFL") : false;
  };
  const modeP = (r) => (gfl(r) ? "PQ" : "PV");
  const modeV = (r) => (gfl(r) ? "VdcQ" : "VdcV");
  const tient = (r) => String(r[3]).startsWith("Vdc");
  const batterie = lignes("T_DCDC").find((r) => ile.has(+r[1]) && String(r[3]).toLowerCase() === "batterie");
  const avant = vscs.map((r) => r[3]);
  // un statisme fixe aussi le niveau de tension ; un formeur offshore
  // (Vtheta) ne le fixe pas, et sa commande ne se change pas
  const statisme = vscs.some((r) => r[3] === "droop");
  if (vscs.length === 1 && batterie) vscs[0][3] = modeP(vscs[0]);
  else if (!vscs.some(tient) && !batterie && !statisme) {
    // aucune reference : un convertisseur, et un seul, la devient. Plusieurs
    // references a des noeuds distincts sont admises (CIGRE B4) : on n'en
    // retire jamais une qu'on a posee
    const ref = vscs.find((r) => r[3] !== "Vtheta");
    if (ref) ref[3] = modeV(ref);
  }
  const tenus = new Set(vscs.filter(tient).map((r) => +r[2]));
  if (!tenus.size && batterie && !statisme) tenus.add(+batterie[1]);
  if (tenus.size) for (const r of lignes("T_BusDC")) if (ile.has(+r[0])) r[1] = tenus.has(+r[0]) ? 1 : 2;
  const ref = vscs.find(tient);
  const change = vscs.filter((r, k) => r[3] !== avant[k]);
  if (!change.length) return null;
  return (
    "commande accordee au reseau continu : " +
    change.map((r) => `convertisseur du bus ${r[1]} en ${r[3]}`).join(", ") +
    (ref ? ` ; la tension continue est tenue par le convertisseur du bus ${ref[1]}` : " ; la batterie tient la tension continue")
  );
}

/** Raccorder le convertisseur d'un bus alternatif a un bus continu.
 *
 *  La ligne de T_VSCs EST la declaration : c'est elle qui fait du GFM3 de ce
 *  bus un GFM3DC, dont le condensateur n'est plus un etat interne mais le bus
 *  continu lui-meme. On aligne au passage la colonne Vdc_V du convertisseur
 *  sur le Un_V du noeud -- l'assembleur exige qu'ils soient egaux, et laisser
 *  l'utilisateur decouvrir ca par un message d'erreur serait avare. */
function poserVSC(busAC, busDC) {
  // Un convertisseur a bus continu IDEAL (GFM1, GFM2, GFL1) n'a pas de port
  // continu : on le passe d'abord dans sa variante a bus continu dynamique,
  // qui est la seule qu'on sache raccorder.
  for (const [de, vers] of [["T_GFM1s", "T_GFM3s"], ["T_GFM2s", "T_GFM4s"], ["T_GFL1s", "T_GFL2s"]]) {
    const i = lignes(de).findIndex((r) => +r[1] === busAC);
    if (i >= 0) {
      const l = lignes(de).splice(i, 1)[0];
      l[0] = nouvelId(vers);
      lignes(vers).push(l);
    }
  }
  const conv = convertisseurDuBus(busAC);
  if (!conv) {
    etat(`bus ${busAC} : aucun GFM ou GFL a raccorder`, true);
    return;
  }
  if (lignes("T_VSCs").some((r) => +r[1] === busAC)) {
    etat(`bus ${busAC} : son convertisseur est deja raccorde`, true);
    return;
  }
  // Refuser tout de suite plutot que de laisser l'assembleur le dire : le
  // premier equipement pose devient la reference d'angle, si bien que le
  // geste naturel -- poser un convertisseur, puis le raccorder -- tombe
  // droit sur cette regle.
  const hote = lignes("Y_DER").find((r) => +r[0] === busAC);
  if (hote && +hote[1] === 1) {
    etat(
      `bus ${busAC} : c'est le slack. Un convertisseur tient exactement deux ` +
      `grandeurs, et sur la reference d'angle les deux puissances sont libres.`,
      true
    );
    return;
  }
  // Le premier raccorde tient la tension continue, les suivants leur
  // puissance : une ile continue a besoin d'exactement une reference.
  const tient = lignes("T_VSCs").some((r) => String(r[3]).startsWith("Vdc"));
  const mode = conv.famille === "GFL2" ? (tient ? "PQ" : "VdcQ") : tient ? "PV" : "VdcV";
  const unDC = (lignes("T_BusDC").find((r) => +r[0] === busDC) || [0, 0, 40e3])[2];
  lignes(conv.table)[conv.index][3] = unDC;
  lignes("T_VSCs").push([nouvelId("T_VSCs"), busAC, busDC, mode]);
  // la regle ci-dessus compte les references de TOUT le reseau ; c'est par
  // ile continue qu'il en faut une, et une batterie peut la fournir
  const msg = accorderIleDC(busDC);
  selectionner("T_VSCs", lignes("T_VSCs").length - 1);
  modifie();
  if (msg) etat(msg);
}

/** Les outils du cote continu : la source qu'ils posent derriere un DC/DC. */
const OUTILS_COTE_DC = { cote_batterie: "batterie", cote_pv: "PV" };

/** Poser une batterie ou un PV du cote continu d'un convertisseur.
 *
 *  ``b`` est un bus continu (``dc``) ou le bus alternatif d'un GFM/GFL. Dans
 *  ce cas, on suit son raccordement (``T_VSCs``) jusqu'a son bus continu ; un
 *  convertisseur a bus continu IDEAL recoit d'abord un bus continu, avec son
 *  condensateur, et passe dans sa variante a port continu (``poserVSC``). La
 *  source est un DC/DC : le PV injecte sa puissance (MPPT), la batterie tient
 *  la tension par statisme. */
function poserSourceCoteDC(b, dc, source) {
  let bDC = dc ? b : null;
  if (!dc) {
    const vsc = lignes("T_VSCs").find((r) => +r[1] === b);
    if (vsc) bDC = +vsc[2];
    else {
      const aConvertisseur = ["T_GFM1s", "T_GFM2s", "T_GFM3s", "T_GFM4s", "T_GFL1s", "T_GFL2s"]
        .some((t) => lignes(t).some((r) => +r[1] === b));
      if (!aConvertisseur) {
        etat(`bus ${b} : cliquer un bus continu ou un GFM/GFL`, true);
        return;
      }
      bDC = nouveauNumero();
      const unDC = lignes("T_BusDC").length ? lignes("T_BusDC")[0][2] : 40e3;
      lignes("T_BusDC").push([bDC, 1, unDC, 1.0, 0]);
      // le bus continu dans le prolongement de la face continue : du cote
      // libre du bus alternatif, a deux longueurs de symbole
      const svg = $("#canevas");
      const q = versPixels(svg, positionBus(b));
      const d = directionLibre(svg, b);
      E.positions[String(bDC)] = placeLibre(versUnite(svg, q[0] + d[0] * 120, q[1] + d[1] * 120), bDC);
      lignes("T_ShuntsDC").push([nouvelId("T_ShuntsDC"), bDC, 0.2]);
      poserVSC(b, bDC);
      if (!lignes("T_VSCs").some((r) => +r[1] === b)) return; // refuse : poserVSC l'a dit
    }
  }
  const Sb = E.tables.Sb;
  lignes("T_DCDC").push([nouvelId("T_DCDC"), bDC, Sb * 0.5, source, Sb * 0.1]);
  const msg = accorderIleDC(bDC);
  selectionner("T_DCDC", lignes("T_DCDC").length - 1);
  modifie();
  etat(`${SOURCES[source].titre} posee au bus continu ${bDC}` + (msg ? ` ; ${msg}` : ""));
}

/** Un convertisseur DC/DC : PV par defaut (puissance constante) ; ses
 *  proprietes le passent en batterie (statisme de tension). */
function poserDCDC(bus) {
  lignes("T_DCDC").push([nouvelId("T_DCDC"), bus, E.tables.Sb * 0.1, "PV", E.tables.Sb * 0.05]);
  const msg = accorderIleDC(bus);
  selectionner("T_DCDC", lignes("T_DCDC").length - 1);
  modifie();
  if (msg) etat(msg);
}

function poserChargeDC(bus) {
  lignes("T_LoadsDC").push([nouvelId("T_LoadsDC"), bus, E.tables.Sb * 0.02, "P"]);
  selectionner("T_LoadsDC", lignes("T_LoadsDC").length - 1);
  modifie();
}

function poserShuntDC(bus) {
  lignes("T_ShuntsDC").push([nouvelId("T_ShuntsDC"), bus, 0.05]);
  selectionner("T_ShuntsDC", lignes("T_ShuntsDC").length - 1);
  modifie();
}

/** Une place libre autour d'un bus, pour y poser un bus terminal.
 *
 *  Le decalage etait FIXE -- en haut a droite de l'hote -- si bien que deux
 *  equipements poses sur le meme bus s'empilaient exactement : une machine
 *  et un GFM au meme noeud n'en montraient qu'un, l'autre dessous. Les pastilles
 *  de couleur le cachaient a moitie ; les symboles le rendent evident. On
 *  essaie donc huit directions autour de l'hote et l'on prend la premiere
 *  assez loin de tout ce qui est deja pose. */
function placeLibre(p, exclu) {
  // R assez grand pour que les symboles de transformateur, poses au milieu de
  // chaque raccordement, ne se touchent pas quand trois unites partagent un hote
  const R = 0.12, marge = 0.06;
  const occupe = Object.entries(E.positions)
    .filter(([k]) => +k !== exclu)
    .map(([, v]) => v);
  const candidats = [];
  for (let k = 0; k < 8; k++) {
    const a = -Math.PI / 4 - (k * Math.PI) / 4;   // en haut a droite d'abord
    candidats.push([
      Math.min(0.98, Math.max(0.02, p[0] + R * Math.cos(a))),
      Math.min(0.98, Math.max(0.02, p[1] + R * Math.sin(a))),
    ]);
  }
  const libre = candidats.find((c) => occupe.every((q) => Math.hypot(q[0] - c[0], q[1] - c[1]) > marge));
  return libre || candidats[0];
}

function poserBranche(de, vers) {
  if (E.outil === "transfo") {
    const u = (b) => (lignes("Y_DER").find((r) => +r[0] === b) || [0, 0, 20e3])[2];
    lignes("T_Transfos").push([nouvelId("T_Transfos"), de, vers, E.tables.Sb, u(de), u(vers), 0.0, 0.05]);
    selectionner("T_Transfos", lignes("T_Transfos").length - 1);
  } else if (E.outil === "ligne_si") {
    lignes("T_Lines_SI").push([nouvelId("T_Lines_SI"), de, vers, 1.0, 0.5, 0.7, 0]);
    selectionner("T_Lines_SI", lignes("T_Lines_SI").length - 1);
  } else {
    const un = (lignes("Y_DER").find((r) => +r[0] === de) || [0, 0, 20e3])[2];
    lignes("T_Lines_pu").push([nouvelId("T_Lines_pu"), de, vers, 1.0, 0.01, 0.1, 0.0, un, E.tables.Sb]);
    selectionner("T_Lines_pu", lignes("T_Lines_pu").length - 1);
  }
  modifie();
}

function supprimerSelection() {
  if (!E.sel) return;
  const { table, index } = E.sel;
  if (table === "__instance") {
    E.reseau.instances.splice(index, 1);
    E.sel = null;
    modifie();
    return;
  }
  if (table === "Y_DER") {
    const b = +lignes(table)[index][0];
    for (const t of TABLES_BRANCHE) E.tables[t] = lignes(t).filter((r) => +r[1] !== b && +r[2] !== b);
    for (const spec of Object.values(E.defauts.equipements))
      E.tables[spec.table] = lignes(spec.table).filter((r) => +r[1] !== b);
    // le raccordement continu part avec le convertisseur qui le portait
    E.tables.T_VSCs = lignes("T_VSCs").filter((r) => +r[1] !== b);
    delete E.positions[String(b)];
  }
  if (table === "T_BusDC") {
    const b = +lignes(table)[index][0];
    E.tables.T_LinesDC = lignes("T_LinesDC").filter((r) => +r[1] !== b && +r[2] !== b);
    E.tables.T_ShuntsDC = lignes("T_ShuntsDC").filter((r) => +r[1] !== b);
    E.tables.T_DCDC = lignes("T_DCDC").filter((r) => +r[1] !== b);
    E.tables.T_LoadsDC = lignes("T_LoadsDC").filter((r) => +r[1] !== b);
    E.tables.T_VSCs = lignes("T_VSCs").filter((r) => +r[2] !== b);
    delete E.positions[String(b)];
  }
  lignes(table).splice(index, 1);
  E.sel = null;
  modifie();
}

// --- inspecteur -------------------------------------------------------
const UNITES = {
  Un_V: ["kV", 1e3], Sn_VA: ["MVA", 1e6], Sb_VA: ["MVA", 1e6], Vdc_V: ["kV", 1e3],
  Pg_W: ["MW", 1e6], Qg_var: ["Mvar", 1e6], Pl_W: ["MW", 1e6], Ql_var: ["Mvar", 1e6],
  P_W: ["MW", 1e6], Q_var: ["Mvar", 1e6], U1_V: ["kV", 1e3], U2_V: ["kV", 1e3], Pn_W: ["MW", 1e6], Pinj_W: ["MW", 1e6],
};

/** Une instance : ce qui la distingue des autres copies de la meme zone. La
 *  zone elle-meme ne s'edite pas ici mais en l'ouvrant, parce qu'une
 *  modification de zone vaut pour TOUTES ses instances. */
function rendreInstance(n, k) {
  const inst = E.reseau && E.reseau.instances[k];
  if (!inst) return;
  n.append(el("div", { class: "mono doux" }, `instance ${k} de la zone ${inst.zone}`));
  const champ = (cle, titre, aide) => {
    n.append(el("label", {}, titre));
    const c = el("input", { type: "number", step: "1", value: inst[cle] });
    c.addEventListener("change", () => {
      const v = parseInt(c.value, 10);
      if (Number.isNaN(v)) return;
      inst[cle] = v;
      modifie();
    });
    n.append(c);
    if (aide) n.append(el("p", { class: "doux" }, aide));
  };
  champ("bus_amont", "Bus amont (transport)");
  champ("prefixe", "Prefixe de numerotation", "Les bus de cette instance sont numerotes a partir de lui.");
  champ("K", "K -- departs identiques", "Puissances multipliees par K, impedances divisees par K.");
  const b = el("button", {}, `Ouvrir la zone ${inst.zone}`);
  b.addEventListener("click", () => ouvrirVue(inst.zone));
  n.append(b);
}

// La topologie d'un cable continu : colonnes optionnelles Config et Circuits
// de T_LinesDC (g2egrid.tables.cable_dc). Les parametres lineiques sont alors
// lus PAR POLE et ramenes au conducteur equivalent par l'assembleur.
function champsTopologieDC(ligne) {
  const configs = (E.defauts && E.defauts.configs_dc) || ["equivalent"];
  const cfg = el("select", {});
  for (const v of configs)
    cfg.append(el("option", { value: v, selected: (ligne[7] || "equivalent") === v ? "" : null }, v));
  const nb = el("input", { type: "number", min: "1", step: "1", value: ligne[8] || 1 });
  const ecrire = () => {
    ligne[7] = cfg.value;
    ligne[8] = Math.max(1, parseInt(nb.value, 10) || 1);
    modifie();
  };
  cfg.addEventListener("change", ecrire);
  nb.addEventListener("change", ecrire);
  return [
    el("label", { title: "equivalent : colonnes du conducteur pole a pole ; sinon : parametres par pole" }, "Config"),
    cfg,
    el("label", {}, "Circuits"),
    nb,
  ];
}

function rendreInspecteur() {
  const n = $("#inspecteur");
  n.textContent = "";
  if (!E.sel) {
    n.append(el("p", { class: "doux" }, "Rien de selectionne."));
    n.append(el("p", { class: "doux" }, "Poser un equipement sur un bus lui cree son bus terminal et son transfo de raccordement."));
    return;
  }
  const { table, index } = E.sel;
  if (table === "__instance") {
    rendreInstance(n, index);
    return;
  }
  const ligne = lignes(table)[index];
  if (!ligne) {
    E.sel = null;
    return;
  }
  n.append(el("div", { class: "mono doux" }, `${table} [${index}]`));

  // le disjoncteur de l'element (pour Y_DER : la charge du bus, s'il en a une)
  if (OUVRABLES.includes(table) && !(table === "Y_DER" && !(+ligne[7] || +ligne[8]))) {
    const id = +ligne[0];
    const cb = el("input", { type: "checkbox", checked: estOuvertId(table, id) ? null : "" });
    cb.addEventListener("change", () => basculerDisjoncteur(table, id, !cb.checked));
    n.append(el("label", { class: "case-service" }, cb,
      table === "Y_DER" ? " charge en service" : " en service (disjoncteur ferme)"));
    if (horsTension(table === "Y_DER" ? id : +ligne[1]))
      n.append(el("p", { class: "doux" }, "hors tension : son ilot n'a plus de formeur de reseau"));
  }

  COLS(table).forEach((col, j) => {
    const lab = el("label", {}, col);
    let champ;
    if (col === "Type" && table === "Y_DER") {
      champ = el("select", {});
      for (const [v, t] of Object.entries(E.defauts.types_bus))
        champ.append(el("option", { value: v, selected: +ligne[j] === +v ? "" : null }, `${v} - ${t}`));
      champ.addEventListener("change", () => {
        ligne[j] = +champ.value;
        modifie();
      });
    } else if (col === "Mode" && table === "T_VSCs") {
      // Les modes possibles dependent de la CLASSE du convertisseur : un
      // grid-forming tient la tension alternative, un grid-following le
      // reactif. Proposer les autres reviendrait a laisser composer un
      // reseau que l'assembleur refusera.
      const conv = convertisseurDuBus(+ligne[1]);
      const modes = (conv && (E.defauts.modes_vsc || {})[conv.famille]) || [];
      champ = el("select", {});
      if (!modes.length) champ.append(el("option", { value: "" }, "(aucun convertisseur sur ce bus)"));
      for (const v of modes)
        champ.append(el("option", { value: v, selected: ligne[j] === v ? "" : null }, v));
      champ.value = ligne[j];
      champ.addEventListener("change", () => {
        ligne[j] = champ.value;
        const msg = accorderModeVSC(ligne);
        modifie();
        if (msg) etat(msg);
      });
    } else if (col === "Type" && table === "T_BusDC") {
      champ = el("select", {});
      for (const [v, t] of Object.entries({ 1: "reference Vdc", 2: "P impose" }))
        champ.append(el("option", { value: v, selected: +ligne[j] === +v ? "" : null }, `${v} - ${t}`));
      champ.addEventListener("change", () => {
        ligne[j] = +champ.value;
        modifie();
      });
    } else if (col === "Type" && (table === "T_DCDC" || table === "T_LoadsDC")) {
      const choix = table === "T_DCDC" ? { PV: "PV (puissance constante)", batterie: "batterie (statisme de tension)" } : { P: "puissance constante", R: "resistance" };
      champ = el("select", {});
      for (const [v, t] of Object.entries(choix))
        champ.append(el("option", { value: v, selected: String(ligne[j]).toLowerCase() === v.toLowerCase() ? "" : null }, t));
      champ.addEventListener("change", () => {
        ligne[j] = champ.value;
        // une batterie peut tenir la tension continue, un PV non : la
        // commande des convertisseurs de l'ile suit
        const msg = table === "T_DCDC" ? accorderIleDC(+ligne[1]) : null;
        modifie();
        if (msg) etat(msg);
      });
    } else if (col === "Modele" && table === "T_Blocs") {
      champ = el("select", {});
      for (const v of Object.keys(E.tables.MODELES_BLOCS || {}))
        champ.append(el("option", { value: v, selected: ligne[j] === v ? "" : null }, v));
      champ.addEventListener("change", () => {
        ligne[j] = champ.value;
        modifie();
      });
      const ouvrir = el("button", { type: "button", title: "Ouvrir le modele dans l'onglet Blocs" }, "editer");
      ouvrir.addEventListener("click", () => {
        B.modele = ligne[j];
        activerOnglet("blocs");
      });
      lab.append(" ", ouvrir);
    } else if (col === "Type" && table === "T_Loads") {
      champ = el("select", {});
      for (const v of ["P", "I"])
        champ.append(el("option", { value: v, selected: ligne[j] === v ? "" : null }, v));
      champ.addEventListener("change", () => {
        ligne[j] = champ.value;
        modifie();
      });
    } else {
      const u = UNITES[col];
      if (u) lab.append(el("span", { class: "doux" }, "  [" + u[0] + "]"));
      champ = el("input", { type: "number", step: "any", value: u ? +ligne[j] / u[1] : ligne[j] });
      champ.addEventListener("change", () => {
        const v = parseFloat(champ.value);
        if (Number.isNaN(v)) return;
        ligne[j] = u ? v * u[1] : v;
        modifie();
      });
    }
    n.append(lab, champ);
  });
  if (table === "T_LinesDC") n.append(...champsTopologieDC(ligne));

  // le statisme Vdc-P d'une station en mode droop : sur sa ligne de T_VSCs
  // comme sur son convertisseur
  const vscDroop =
    table === "T_VSCs" ? ligne : lignes("T_VSCs").find((r) => /^T_GF[ML]\ds$/.test(table) && +r[1] === +ligne[1]);
  if (/^T_GF[ML]\ds$/.test(table)) n.append(blocCommande(table, index));
  if (vscDroop && /^droop/.test(String(vscDroop[3]))) n.append(blocStatisme(vscDroop));

  // variante, dynamiques et parametres de l'equipement (SM, GFM, GFL)
  const cle = familleReglages(table);
  if (table === "T_GFL2s" && (E.tables.T_VSCs || []).some((r) => +r[1] === +ligne[1]))
    n.append(blocModeleStation(+ligne[1]));
  if (cle) {
    const bus = +ligne[1];
    n.append(blocDynamiques(cle, table, index, bus));
    n.append(blocParametres(cle, bus, table));
  }
  if (classesModele(table, ligne).length)
    n.append(el("button", { class: "btn-fiche", title: "ou double-clic sur l'element du schema",
      onclick: () => ouvrirFiche(table, index) }, "Fiche : schemas et equations"));
  if (typeof TABLES_EN_BLOCS !== "undefined" && TABLES_EN_BLOCS.test(table)) {
    const sortie = el("div", {});
    n.append(el("div", { class: "bloc-en-blocs" }, ...boutonsEnBlocs(table, index, sortie)), sortie);
  }
  if (typeof sectionEquations === "function" && TABLES_EQUATIONS.test(table)) n.append(sectionEquations(table, index, "inspecteur"));
  for (const c of classesModele(table, ligne)) n.append(blocEquations(c));
}

// --- commande d'un convertisseur : ce qu'il tient, de chaque cote -----------
/** Les choix du cote continu et du cote alternatif, et ce qu'ils deviennent
 *  dans les tables : le mode du VSC (``T_VSCs``), le type des bus (PV : tient
 *  V, PQ : tient Q, slack : forme son ile), ``statisme_qv`` (statisme Q-V) et
 *  ``kdroop`` (statismes continus) dans ``PARAMS_PAR_BUS``. */
const COMMANDE_DC = {
  Vdc: "tension continue Vdc",
  P: "puissance P",
  droop: "statisme P-Vdc",
  droopI: "statisme Idc-Vdc",
};
const COMMANDE_AC = {
  V: "tension alternative V (bus PV)",
  Q: "puissance reactive Q (bus PQ)",
  QV: "statisme Q-V (bus PQ)",
  Vtheta: "forme le reseau : V et angle (slack)",
};
const LOIS_GFM_UI = {
  statisme: "statisme P-f (1er ordre)",
  vsm: "machine synchrone virtuelle (statisme filtre, 2e ordre)",
  dvoc: "oscillateur virtuel (dVOC)",
  matching: "appariement w = Vdc (matching)",
};

/** Le mode de T_VSCs qui realise le couple (cote continu, cote alternatif). */
function modeCommande(dc, ac) {
  if (ac === "Vtheta") return "Vtheta";
  const tientV = ac === "V";
  if (dc === "Vdc") return tientV ? "VdcV" : "VdcQ";
  if (dc === "P") return tientV ? "PV" : "PQ";
  return dc; // droop, droopI : le cote alternatif est le type du bus
}

function commandeActuelle(busAC) {
  const vsc = lignes("T_VSCs").find((r) => +r[1] === busAC);
  const ac = lignes("Y_DER").find((r) => +r[0] === busAC);
  const reg = (E.tables.PARAMS_PAR_BUS || {})[busAC] || {};
  const mode = vsc ? String(vsc[3]) : null;
  let cAC = ac && +ac[1] === 1 ? "Vtheta" : ac && +ac[1] === 2 ? "V" : +reg.statisme_qv > 0 ? "QV" : "Q";
  let cDC = !mode ? null : mode.startsWith("Vdc") ? "Vdc" : mode === "droop" || mode === "droopI" ? mode : mode === "Vtheta" ? null : "P";
  return { vsc, ac, reg, cDC, cAC };
}

function blocCommande(table, index) {
  const ligne = lignes(table)[index];
  const busAC = +ligne[1];
  const gfm = table.startsWith("T_GFM");
  const d = el("details", { open: "", class: "bloc-commande" });
  d.append(el("summary", {}, "Commande : ce que tient le convertisseur"));
  const { vsc, cDC, cAC } = commandeActuelle(busAC);
  const conv = convertisseurDuBus(busAC);
  const admis = vsc && conv ? (E.defauts.modes_vsc || {})[conv.famille] || [] : null;

  const choix = (titre, options, valeur, possible, ecrire, aide) => {
    const sel = el("select", { title: aide || "" });
    for (const [v, t] of Object.entries(options)) {
      const ok = possible(v);
      sel.append(el("option", { value: v, selected: v === valeur ? "" : null, disabled: ok ? null : "" },
        ok ? t : `${t} (impossible ici)`));
    }
    sel.value = valeur;
    sel.addEventListener("change", () => ecrire(sel.value));
    d.append(el("label", {}, titre), sel);
    return sel;
  };

  // cote continu
  if (vsc && cAC === "Vtheta") {
    d.append(el("label", {}, "cote continu"),
      el("p", { class: "doux" }, "puissance libre : celle que demande l'ile alternative qu'il forme ; sa tension continue est tenue ailleurs dans l'ile continue"));
  } else if (vsc) {
    choix("cote continu", COMMANDE_DC, cDC,
      (v) => cAC !== "Vtheta" && admis.includes(modeCommande(v, cAC === "V" ? "V" : "Q")),
      (v) => appliquerCommande(busAC, v, cAC),
      "ce que le convertisseur regle sur son noeud continu ; un seul convertisseur (ou une batterie) tient la tension d'une ile continue");
  } else {
    d.append(el("label", {}, "cote continu"),
      el("p", { class: "doux" }, "bus continu ideal (source sans dynamique). Pour le detailler, poser une batterie ou un PV de son cote continu (outils de la palette)."));
  }
  // cote alternatif
  choix("cote alternatif", COMMANDE_AC, cAC,
    (v) => {
      if (v === "QV") return !!vsc && cDC !== null;
      if (!vsc) return v !== "QV";
      if (v === "Vtheta") return admis.includes("Vtheta");
      return admis.includes(modeCommande(cDC || "P", v === "V" ? "V" : "Q")) || ["droop", "droopI"].includes(cDC);
    },
    (v) => appliquerCommande(busAC, cDC, v),
    "ce que le convertisseur regle sur son bus alternatif : une tension (PV), un reactif (PQ), un statisme Q-V, ou tout le reseau qu'il forme (slack)");
  if (cAC === "QV") {
    const reg = E.tables.PARAMS_PAR_BUS[busAC];
    const inp = el("input", { type: "number", step: "any", value: +(100 * reg.statisme_qv).toPrecision(6) });
    inp.addEventListener("change", () => {
      const v = parseFloat(inp.value);
      if (!(v > 0)) return;
      reg.statisme_qv = v / 100;
      modifie();
    });
    d.append(el("label", { title: "pu de tension par pu de reactif, base Sn ; consigne : V du bus, reactif a V = consigne : Qg du bus" }, "nq [%]"), inp);
  }
  // loi de synchronisation d'un grid-forming
  if (gfm) {
    const reg = ((E.tables.PARAMS_PAR_BUS || {})[busAC]) || {};
    const loiReg = String(reg.loi || "statisme").toLowerCase();
    const loi = loiReg !== "statisme" ? loiReg : table === "T_GFM2s" || table === "T_GFM4s" ? "vsm" : "statisme";
    const dcDyn = table === "T_GFM3s" || table === "T_GFM4s";
    choix("loi de synchronisation", LOIS_GFM_UI, loi,
      (v) => (v === "dvoc" ? !dcDyn : v === "matching" ? dcDyn : true),
      (v) => appliquerLoiGFM(table, index, v),
      "comment le grid-forming fixe sa frequence : statisme P-f, statisme filtre (inertie virtuelle, VSM), oscillateur virtuel, ou appariement sur la tension continue");
  }
  return d;
}

/** Ecrire le couple (cote continu, cote alternatif) dans les tables. */
function appliquerCommande(busAC, dc, ac) {
  const { vsc, ac: bus } = commandeActuelle(busAC);
  E.tables.PARAMS_PAR_BUS = E.tables.PARAMS_PAR_BUS || {};
  const reg = (E.tables.PARAMS_PAR_BUS[busAC] = E.tables.PARAMS_PAR_BUS[busAC] || {});
  const msgs = [];
  if (vsc) {
    vsc[3] = modeCommande(dc || "P", ac);
    const m = accorderModeVSC(vsc);
    if (m) msgs.push(m);
  }
  if (bus && ac !== "Vtheta") {
    const type = ac === "V" ? 2 : 3;
    if (+bus[1] !== type) {
      bus[1] = type;
      msgs.push(`bus ${busAC} en ${type === 2 ? "PV" : "PQ"}`);
    }
  } else if (bus && ac === "Vtheta" && +bus[1] !== 1) {
    bus[1] = 1;
    msgs.push(`bus ${busAC} en slack`);
  }
  if (ac === "QV") {
    if (!(+reg.statisme_qv > 0)) {
      reg.statisme_qv = 0.05;
      msgs.push("statisme Q-V nq = 5 %");
    }
  } else delete reg.statisme_qv;
  modifie();
  rendreInspecteur();
  if (msgs.length) etat(msgs.join(" ; "));
}

/** La loi de synchronisation d'un GFM : la table (1er ou 2e ordre) et le
 *  reglage ``loi`` (dVOC sur bus continu ideal, matching sur bus dynamique). */
function appliquerLoiGFM(table, index, loi) {
  const ligne = lignes(table)[index];
  const busAC = +ligne[1];
  E.tables.PARAMS_PAR_BUS = E.tables.PARAMS_PAR_BUS || {};
  const reg = (E.tables.PARAMS_PAR_BUS[busAC] = E.tables.PARAMS_PAR_BUS[busAC] || {});
  const dcDyn = table === "T_GFM3s" || table === "T_GFM4s";
  const second = loi === "vsm";
  if (loi === "dvoc" || loi === "matching") reg.loi = loi;
  else delete reg.loi;
  // dVOC et matching partent du 1er ordre ; la VSM est le 2e ordre
  const cible = dcDyn ? (second ? "T_GFM4s" : "T_GFM3s") : second ? "T_GFM2s" : "T_GFM1s";
  if (cible !== table) changerTable(table, index, cible);
  else {
    modifie();
    rendreInspecteur();
  }
  etat(`GFM du bus ${busAC} : ${LOIS_GFM_UI[loi]}`);
}

// --- statisme Vdc-P ---------------------------------------------------------
/** Ce qu'un mode de station demande aux tables, applique au moment ou on le
 *  choisit plutot que decouvert au load flow :
 *
 *  * ``droop`` : le noeud continu est de type 2 (sa V_pu est la consigne) et
 *    la station recoit un statisme de 5 % s'il n'en a pas ;
 *  * ``Vdc...`` : la station tient son noeud continu, de type 1 ;
 *  * ``Vtheta`` : la station forme son ile, son bus alternatif est le slack ;
 *    en quittant ce mode, il redevient PV (GFM) ou PQ (GFL). */
function accorderModeVSC(vsc) {
  const mode = String(vsc[3]);
  const busAC = +vsc[1], busDC = +vsc[2];
  const dc = lignes("T_BusDC").find((r) => +r[0] === busDC);
  const ac = lignes("Y_DER").find((r) => +r[0] === busAC);
  const conv = convertisseurDuBus(busAC);
  const msgs = [];
  if (mode === "droop" || mode === "droopI") {
    if (dc && +dc[1] !== 2) {
      dc[1] = 2;
      msgs.push(`noeud continu ${busDC} en type 2 (sa tension ${fmt(+dc[3], 4)} pu devient la consigne du statisme)`);
    }
    E.tables.PARAMS_PAR_BUS = E.tables.PARAMS_PAR_BUS || {};
    const reg = (E.tables.PARAMS_PAR_BUS[busAC] = E.tables.PARAMS_PAR_BUS[busAC] || {});
    if (!(+reg.kdroop > 0)) {
      reg.kdroop = 0.05;
      msgs.push("statisme k = 5 % (pu de tension continue par pu de puissance)");
    }
  } else if (mode.startsWith("Vdc")) {
    if (dc && +dc[1] !== 1) {
      dc[1] = 1;
      msgs.push(`noeud continu ${busDC} en type 1 (tenu par la station)`);
    }
  }
  if (ac) {
    if (mode === "Vtheta" && +ac[1] !== 1) {
      ac[1] = 1;
      msgs.push(`bus ${busAC} en slack : la station forme son ile`);
    } else if (mode !== "Vtheta" && +ac[1] === 1) {
      ac[1] = conv && conv.famille.startsWith("GFL") ? 3 : 2;
      msgs.push(`bus ${busAC} en ${+ac[1] === 2 ? "PV" : "PQ"}`);
    }
  }
  return msgs.length ? `mode ${mode} : ` + msgs.join(" ; ") : null;
}

/** Le statisme d'une station : k, P_ref, Vdc_ref, et sa caracteristique
 *  P(Vdc), avec le point de fonctionnement du dernier load flow. */
function blocStatisme(vsc) {
  const busAC = +vsc[1], busDC = +vsc[2];
  const ac = lignes("Y_DER").find((r) => +r[0] === busAC);
  const dc = lignes("T_BusDC").find((r) => +r[0] === busDC);
  const conv = convertisseurDuBus(busAC);
  const sn = conv ? +lignes(conv.table)[conv.index][2] : E.tables.Sb;
  E.tables.PARAMS_PAR_BUS = E.tables.PARAMS_PAR_BUS || {};
  const reg = () => (E.tables.PARAMS_PAR_BUS[busAC] = E.tables.PARAMS_PAR_BUS[busAC] || {});
  const d = el("details", { open: "", class: "bloc-statisme" });
  const courant = String(vsc[3]) === "droopI";
  d.append(el("summary", {}, courant ? "Statisme Idc-Vdc" : "Statisme Vdc-P"));
  const champ = (texte, valeur, ecrire, aide) => {
    const lab = el("label", { title: aide }, texte);
    const inp = el("input", { type: "number", step: "any", value: +(+valeur).toPrecision(6) });
    inp.addEventListener("change", () => {
      const v = parseFloat(inp.value);
      if (Number.isNaN(v)) return;
      ecrire(v);
      modifie();
    });
    d.append(lab, inp);
  };
  champ("k [%]", 100 * (+reg().kdroop || 0), (v) => { reg().kdroop = v / 100; },
    "pu de tension continue par pu de puissance, sur la puissance Sn de la station");
  champ("P_ref [MW]", ac ? +ac[5] / 1e6 : 0, (v) => { if (ac) ac[5] = v * 1e6; },
    "puissance injectee dans l'alternatif a Vdc = Vdc_ref (Pg du bus de la station)");
  champ("Vdc_ref [pu]", dc ? +dc[3] : 1, (v) => { if (dc) dc[3] = v; },
    "consigne de tension continue (V_pu du noeud continu)");
  d.append(caracteristiqueStatisme(busAC, busDC, sn));
  d.append(el("p", { class: "doux" },
    (courant
      ? "Idc = P_ref/Vdc_ref + (Sn/k).(Vdc - Vdc_ref), soit P = Vdc.Idc : la station injecte davantage quand la tension continue monte. "
      : "P = P_ref + (Sn/k).(Vdc - Vdc_ref) : la station injecte davantage quand la tension continue monte. ") +
    "Le balayage des statismes est dans l'onglet Limites et impedances."));
  return d;
}

function caracteristiqueStatisme(busAC, busDC, sn) {
  const ac = lignes("Y_DER").find((r) => +r[0] === busAC);
  const dc = lignes("T_BusDC").find((r) => +r[0] === busDC);
  const k = +(((E.tables.PARAMS_PAR_BUS || {})[busAC] || {}).kdroop || 0);
  const pref = ac ? +ac[5] / 1e6 : 0, vref = dc ? +dc[3] : 1;
  const W = 260, H = 150, ML = 46, MB = 26, MT = 8, MR = 8;
  const s = el("svg", { viewBox: `0 0 ${W} ${H}`, width: "100%", class: "carac-statisme" });
  if (!(k > 0)) return s;
  const pente = sn / 1e6 / k; // MW par pu
  const v0 = vref - 0.05, v1 = vref + 0.05;
  const p = (v) => pref + pente * (v - vref);
  let pts = [[v0, p(v0)], [v1, p(v1)]];
  const lf = E.res.lf;
  const vdc = lf && lf.bus_dc ? (lf.bus_dc.find((x) => x.id === busDC) || {}).V : undefined;
  const pac = lf ? (lf.bus.find((x) => x.id === busAC) || {}).P : undefined;
  const op = vdc !== undefined && pac !== undefined ? [vdc, pac * E.tables.Sb / 1e6] : null;
  const ys = pts.map((q) => q[1]).concat(op ? [op[1]] : []);
  const ymin = Math.min(...ys), ymax = Math.max(...ys), dy = ymax - ymin || 1;
  const X = (v) => ML + ((v - v0) / (v1 - v0)) * (W - ML - MR);
  const Y = (y) => MT + (1 - (y - ymin) / dy) * (H - MT - MB);
  s.append(el("line", { x1: ML, y1: H - MB, x2: W - MR, y2: H - MB, stroke: "#98a2b3" }));
  s.append(el("line", { x1: ML, y1: MT, x2: ML, y2: H - MB, stroke: "#98a2b3" }));
  s.append(el("polyline", { points: pts.map(([v, y]) => `${X(v)},${Y(y)}`).join(" "), fill: "none", stroke: "#7a4fa3", "stroke-width": 2 }));
  s.append(el("circle", { cx: X(vref), cy: Y(pref), r: 3, fill: "#7a4fa3" }));
  if (op) {
    const c = el("circle", { cx: X(op[0]), cy: Y(op[1]), r: 4.5, fill: "none", stroke: "#bf2617", "stroke-width": 2 });
    infobulle(c, `load flow : Vdc = ${fmt(op[0], 4)} pu, P = ${fmt(op[1], 4)} MW`);
    s.append(c);
  }
  s.append(el("text", { x: (ML + W - MR) / 2, y: H - 6, "text-anchor": "middle", class: "symbole-indice", fill: "#344054" }, "Vdc [pu]"));
  s.append(el("text", { x: ML - 4, y: Y(ymax) + 4, "text-anchor": "end", class: "symbole-indice", fill: "#344054" }, fmt(ymax, 3)));
  s.append(el("text", { x: ML - 4, y: Y(ymin), "text-anchor": "end", class: "symbole-indice", fill: "#344054" }, fmt(ymin, 3)));
  s.append(el("text", { x: X(v0), y: H - MB + 12, "text-anchor": "start", class: "symbole-indice", fill: "#344054" }, fmt(v0, 3)));
  s.append(el("text", { x: X(v1), y: H - MB + 12, "text-anchor": "end", class: "symbole-indice", fill: "#344054" }, fmt(v1, 3)));
  return s;
}

// --- equations du modele ----------------------------------------------------
// La classe que l'assembleur posera pour cette ligne -- la meme regle que
// g2egrid.tables.assembler : variante a transfo integre d'une machine, variante
// a port continu d'un convertisseur declare dans T_VSCs...
function classesModele(table, ligne) {
  const vsc = (E.tables.T_VSCs || []).some((r) => +r[1] === +ligne[1]);
  const plein = (v) => v !== undefined && v !== null && v !== "";
  switch (table) {
    case "Y_DER": return +ligne[7] || +ligne[8] ? ["Bus", "RLLoad"] : ["Bus"];
    case "T_Lines_pu": case "T_Lines_SI": return ["Line"];
    case "T_Transfos": return ["Transformer"];
    case "T_SMs": return [plein(ligne[3]) ? "SynchronousMachineTR" : "SynchronousMachine"];
    case "T_GFM1s": return ["GFM1"];
    case "T_GFM2s": return ["GFM2"];
    case "T_GFM3s": return [vsc ? "GFM3DC" : "GFM3"];
    case "T_GFM4s": return [vsc ? "GFM4DC" : "GFM4"];
    case "T_GFL1s": return ["GFL1"];
    case "T_GFL2s": return [vsc ? ({ mmc: "MMC", mmc_energie: "MMCEnergie" }[modeleStation(+ligne[1])] || "GFL2DC") : "GFL2"];
    case "T_IBs": return ["IB"];
    case "T_Loads": return [String(ligne[4]).toUpperCase() === "P" ? "PLoad" : "ILoad"];
    case "T_Shunts": return +ligne[2] < 0 ? ["ReactanceShunt"] : [];
    case "T_BusDC": return ["BusDC"];
    case "T_DCDC": return ["ConvertisseurDCDC"];
    case "T_LoadsDC": return [String(ligne[3]).toUpperCase() === "R" ? "ChargeDC_R" : "ChargeDC_P"];
    case "T_LinesDC": return ["LineDC"];
    default: return [];
  }
}

/** Les equations du modele, telles que SymPy les compile : ce qui s'affiche
 *  est ce qui se calcule. Chargees a l'ouverture seulement. */
function blocEquations(classe) {
  const d = el("details", { class: "avance" }, el("summary", {}, `Equations du modele ${classe}`));
  d.addEventListener("toggle", async () => {
    if (!d.open || d.dataset.charge) return;
    d.dataset.charge = "1";
    try {
      const r = await api("/api/equations/" + classe);
      const bloc = (titre, lignes2, gauche) => {
        if (!lignes2.length) return;
        d.append(el("div", { class: "doux", style: "font-size:12px;margin-top:6px" }, titre));
        const pre = el("pre", { class: "equations" });
        pre.textContent = lignes2.map((l) => `${gauche(l.variable)} = ${l.texte}`).join("\n\n");
        d.append(pre);
      };
      if (r.resume) d.append(el("p", { class: "doux", style: "font-size:12px" }, r.resume));
      bloc("dx/dt = f(x, z, u)", r.f, (v) => `d${v}/dt`);
      bloc("0 = g(x, z, u)", r.g, () => "0");
      bloc("y = h(x, z, u)", r.h, (v) => v);
      const latex = [...r.f, ...r.g, ...r.h].map((l) => l.latex).join("\n");
      const b = el("button", {}, "Copier en LaTeX");
      b.addEventListener("click", () => navigator.clipboard && navigator.clipboard.writeText(latex));
      const vers = el("button", {}, "Voir dans la documentation");
      vers.addEventListener("click", () => ouvrirDoc("equations_modeles.html#" + classe.toLowerCase()));
      d.append(b, vers, el("p", { class: "doux", style: "font-size:12px" },
        "Parametres : " + r.parametres.join(", ") + "."));
    } catch (err) {
      d.append(el("div", { class: "pb erreur" }, err.message));
    }
  });
  return d;
}

// --- reglage par intention ------------------------------------------------
// Les reglages EFFECTIFS d'une unite : defauts de la famille, puis reglages
// de famille des tables, puis ceux de son bus -- l'ordre de l'assembleur.
function reglagesEffectifs(cle, bus, table) {
  const fam = { SM: "SM_PARAMS", GFM: "GFM_PARAMS", GFL: "GFL_PARAMS" }[cle];
  return Object.assign({}, defautsReglages(cle, bus, table), E.tables[fam] || {},
    (E.tables.PARAMS_PAR_BUS || {})[String(bus)] || {});
}

/** Le modele d'une station GFL a port continu (g2egrid.tables.MODELES_STATION) :
 *  celui du bus, a defaut celui de la famille, a defaut "vsc". */
function modeleStation(bus) {
  const propre = ((E.tables.PARAMS_PAR_BUS || {})[String(bus)] || {}).modele;
  return String(propre || (E.tables.GFL_PARAMS || {}).modele || "vsc").toLowerCase();
}

/** Les reglages par defaut d'un equipement : ceux du MMC pour une station
 *  MMC, ceux de sa famille sinon -- avec, pour un GFM dont la tension
 *  continue est une variable (``table`` T_GFM3s ou T_GFM4s), ceux de
 *  g2egrid.elements.limiteurs.REGLAGES_BUS_CONTINU. */
function defautsReglages(cle, bus, table) {
  const mmc = cle === "GFL" && modeleStation(bus).startsWith("mmc");
  if (mmc) return E.defauts.reglages.MMC || {};
  const dc = E.defauts.reglages_bus_continu || { tables: [], valeurs: {} };
  return Object.assign({}, E.defauts.reglages[cle] || {},
    dc.tables.includes(table) ? dc.valeurs : {});
}

/** Le modele d'une station a port continu : VSC a deux niveaux (filtre LCL)
 *  ou MMC moyen, energie libre ou regulee. */
function blocModeleStation(bus) {
  const noms = { vsc: "VSC deux niveaux (filtre LCL)", mmc: "MMC moyen, energie libre", mmc_energie: "MMC moyen, energie regulee" };
  const choix = el("select", {});
  for (const v of E.defauts.modeles_station || Object.keys(noms))
    choix.append(el("option", { value: v, selected: modeleStation(bus) === v ? "" : null }, noms[v] || v));
  choix.addEventListener("change", () => ecrireParBus(bus, { modele: choix.value }));
  return el("div", { class: "bloc-modele-station" },
    el("label", { title: "g2egrid.tables.MODELES_STATION" }, "Modele de station"), choix);
}

/** La famille de reglages d'une table d'unite, ou null. */
function familleReglages(table) {
  if (table === "T_SMs") return "SM";
  if (/^T_GFM\ds$/.test(table)) return "GFM";
  if (/^T_GFL\ds$/.test(table)) return "GFL";
  return null;
}

// Les variantes d'un GFM sont deux dynamiques : [statisme du 2e ordre, bus continu]
const TABLE_GFM = { "00": "T_GFM1s", "10": "T_GFM2s", "01": "T_GFM3s", "11": "T_GFM4s" };

/** Deplacer une unite d'une table a une autre (GFM1 -> GFM3...) : meme ligne,
 *  nouvel identifiant dans sa nouvelle table, selection conservee. */
function changerTable(table, index, cible) {
  if (table === cible) return;
  const ligne = lignes(table).splice(index, 1)[0];
  ligne[0] = nouvelId(cible);
  lignes(cible).push(ligne);
  selectionner(cible, lignes(cible).length - 1);
  modifie();
}

/** Les dynamiques de l'unite, en cases a cocher.
 *
 *  Une case DECOCHEE supprime sa dynamique : l'etat devient statique
 *  (dx/dt = 0) ou fige (x = x0) -- cle ``dynamiques`` des reglages du bus,
 *  voir Element.reduire. Les variantes qui ont leur propre table (statisme du
 *  2e ordre et bus continu d'un GFM, commande Vdc/Q d'un GFL) deplacent
 *  l'unite de table. */
function blocDynamiques(cle, table, index, bus) {
  const d = el("details", { class: "avance", open: "" }, el("summary", {}, "Variante et dynamiques"));
  const pb = (E.tables.PARAMS_PAR_BUS = E.tables.PARAMS_PAR_BUS || {});
  const dyn = Object.assign({}, ((E.tables[{ SM: "SM_PARAMS", GFM: "GFM_PARAMS", GFL: "GFL_PARAMS" }[cle]] || {}).dynamiques) || {},
    ((pb[String(bus)] || {}).dynamiques) || {});
  const raccorde = lignes("T_VSCs").find((r) => +r[1] === bus);
  const numero = (table.match(/\d/) || ["1"])[0];
  for (const g of (E.defauts.dynamiques || {})[cle] || []) {
    let actif, bloque = false, note = "";
    if (g.table && cle === "GFM") {
      actif = g.cle === "statisme2" ? ["2", "4"].includes(numero) : ["3", "4"].includes(numero);
      if (g.cle === "dc" && raccorde) {
        bloque = true;
        note = ` (raccorde au bus continu ${raccorde[2]})`;
      }
    } else if (g.table && cle === "GFL") {
      actif = numero === "2";
      if (raccorde) {
        bloque = true;
        note = ` (raccorde au bus continu ${raccorde[2]})`;
      }
    } else {
      actif = !g.etats.some((e) => dyn[e] && dyn[e] !== "dynamique");
    }
    const cb = el("input", { type: "checkbox", checked: actif ? "" : null, disabled: bloque ? "" : null });
    cb.addEventListener("change", () => {
      if (g.table && cle === "GFM") {
        const s2 = g.cle === "statisme2" ? cb.checked : ["2", "4"].includes(numero);
        const dc = g.cle === "dc" ? cb.checked : ["3", "4"].includes(numero);
        changerTable(table, index, TABLE_GFM[`${+s2}${+dc}`]);
        return;
      }
      if (g.table && cle === "GFL") {
        changerTable(table, index, cb.checked ? "T_GFL2s" : "T_GFL1s");
        return;
      }
      const reg = Object.assign({}, pb[String(bus)] || {});
      const dd = Object.assign({}, reg.dynamiques || {});
      for (const e of g.etats) {
        if (cb.checked) dd[e] = "dynamique";
        else dd[e] = g.sans;
      }
      // "dynamique" n'a de sens que pour surcharger la famille : sinon on l'omet
      const fam = ((E.tables[{ SM: "SM_PARAMS", GFM: "GFM_PARAMS", GFL: "GFL_PARAMS" }[cle]] || {}).dynamiques) || {};
      for (const e of g.etats) if (dd[e] === "dynamique" && !fam[e]) delete dd[e];
      if (Object.keys(dd).length) reg.dynamiques = dd;
      else delete reg.dynamiques;
      if (Object.keys(reg).length) pb[String(bus)] = reg;
      else delete pb[String(bus)];
      modifie();
    });
    const lab = el("label", { class: "case-service", title: g.etats ? `etats : ${g.etats.join(", ")} ; decoche : ${g.sans === "fige" ? "fige (x = x0)" : "statique (dx/dt = 0)"}` : "" }, cb, " " + g.titre + note);
    d.append(lab);
  }
  d.append(el("p", { class: "doux", style: "font-size:12px;margin:6px 0" },
    `Modele : ${classesModele(table, lignes(table)[index] || []).join(", ")}. Une dynamique decochee devient statique (dx/dt = 0) ou figee (x = x0).`));
  return d;
}

/** Tous les parametres de l'unite, en cellules a remplir : valeurs effectives
 *  (defauts de la classe, puis reglages de famille, puis ceux du bus). Une
 *  cellule modifiee s'ecrit dans PARAMS_PAR_BUS ; les temps de reponse des
 *  boucles et, pour une machine, sa fiche constructeur sont dans le meme bloc. */
function blocParametres(cle, bus, table) {
  const d = el("details", { class: "avance" }, el("summary", {}, "Parametres"));
  const eff = reglagesEffectifs(cle, bus, table);
  const defaut = defautsReglages(cle, bus, table);
  const propres = (E.tables.PARAMS_PAR_BUS || {})[String(bus)] || {};
  const grille = el("div", { class: "grille-fiche" });
  const boucle = (k) => /^(tau|zeta)_/.test(k) && cle !== "SM";
  for (const k of Object.keys(defaut)) {
    if (boucle(k) || typeof defaut[k] !== "number") continue;
    const v = eff[k];
    const inp = el("input", { type: "number", step: "any", value: v === null || v === undefined ? "" : +(+v).toPrecision(6), style: "width:100px" });
    if (k in propres) inp.style.fontWeight = "600";
    inp.title = `defaut : ${defaut[k]}`;
    inp.addEventListener("change", () => {
      const x = parseFloat(inp.value);
      if (!Number.isNaN(x)) ecrireParBus(bus, { [k]: x });
    });
    grille.append(el("label", { title: `defaut : ${defaut[k]}` }, k), inp);
  }
  d.append(grille);
  const raz = el("button", {}, "Revenir aux valeurs par defaut");
  raz.addEventListener("click", () => {
    const pb = E.tables.PARAMS_PAR_BUS || {};
    // les dynamiques, le modele de station et la ligne d'arbre ne sont pas des parametres
    const propre = pb[String(bus)] || {};
    const garde = {};
    for (const k of ["dynamiques", "modele", "arbre"]) if (propre[k] !== undefined) garde[k] = propre[k];
    if (Object.keys(garde).length) pb[String(bus)] = garde;
    else delete pb[String(bus)];
    modifie();
  });
  d.append(raz, el("p", { class: "doux", style: "font-size:12px" },
    "En gras : valeurs propres a cet equipement ; les autres viennent des defauts ou des reglages de famille."));
  if (cle === "GFM" || (cle === "GFL" && !modeleStation(bus).startsWith("mmc"))) d.append(blocReglageBoucles(cle, bus));
  if (cle === "SM") d.append(blocDonneesConstructeur(bus));
  return d;
}

function ecrireParBus(bus, valeurs) {
  const pb = (E.tables.PARAMS_PAR_BUS = E.tables.PARAMS_PAR_BUS || {});
  pb[String(bus)] = Object.assign({}, pb[String(bus)] || {}, valeurs);
  modifie();
}

/** Chaque boucle d'un convertisseur par son temps de reponse et son
 *  amortissement ; les gains qui en sortent sont lus sur le modele. */
function blocReglageBoucles(cle, bus) {
  const d = el("details", { class: "avance", open: "" }, el("summary", {}, "Reglage des boucles (temps de reponse, amortissement)"));
  const t = el("table");
  t.append(el("thead", {}, el("tr", {}, el("th", {}, "boucle"), el("th", {}, "tau [ms]"), el("th", {}, "zeta"), el("th", {}, "gains"))));
  const tb = el("tbody");
  t.append(tb);
  d.append(t, el("p", { class: "doux", style: "font-size:12px" },
    "wn = 3/(zeta.tau) ; les gains se recalculent au placement de poles du modele."));
  api("/api/reglage/boucles", { famille: cle, reglages: reglagesEffectifs(cle, bus), f0: E.tables.f0 || 50 })
    .then((r) => {
      for (const b of r.boucles) {
        const tau = el("input", { type: "number", step: "any", value: +(b.tau * 1e3).toPrecision(4), style: "width:70px" });
        tau.addEventListener("change", () => ecrireParBus(bus, { [b.cle_tau]: parseFloat(tau.value) * 1e-3 }));
        const zeta = b.cle_zeta
          ? el("input", { type: "number", step: "0.05", value: b.zeta, style: "width:60px" })
          : el("span", { class: "doux" }, "1er ordre");
        if (b.cle_zeta) zeta.addEventListener("change", () => ecrireParBus(bus, { [b.cle_zeta]: parseFloat(zeta.value) }));
        const gains = Object.entries(b.gains).map(([k, v]) => `${k} = ${fmt(v, 4)}`).join(", ");
        tb.append(el("tr", {}, el("td", {}, b.nom), el("td", {}, tau), el("td", {}, zeta), el("td", { class: "mono" }, gains)));
      }
    })
    .catch((err) => d.append(el("div", { class: "pb erreur" }, err.message)));
  return d;
}

/** Saisir une machine par sa fiche : reactances et constantes de temps. */
function blocDonneesConstructeur(bus) {
  const d = el("details", { class: "avance" }, el("summary", {}, "Donnees constructeur (Xd, X'd, T'd0...)"));
  const noms = {
    Xd: "Xd", Xq: "Xq", Xl: "Xl", Xd_p: "X'd", Xd_pp: "X''d", Xq_p: "X'q", Xq_pp: "X''q",
    Td0_p: "T'd0 [s]", Td0_pp: "T''d0 [s]", Tq0_p: "T'q0 [s]", Tq0_pp: "T''q0 [s]", Ra: "Ra",
  };
  const champs = {};
  const grille = el("div", { class: "grille-fiche" });
  d.append(grille);
  const f0 = E.tables.f0 || 50;
  api("/api/machine/standard", { valeurs: reglagesEffectifs("SM", bus), f0, sens: "vers_standard" })
    .then((r) => {
      for (const [k, lab] of Object.entries(noms)) {
        champs[k] = el("input", { type: "number", step: "any", value: +r.standard[k].toPrecision(5), style: "width:80px" });
        grille.append(el("label", {}, lab), champs[k]);
      }
      const b = el("button", {}, "Appliquer a cette machine");
      b.addEventListener("click", async () => {
        const valeurs = Object.fromEntries(Object.entries(champs).map(([k, n]) => [k, parseFloat(n.value)]));
        try {
          const p = await api("/api/machine/standard", { valeurs, f0, sens: "depuis_standard" });
          ecrireParBus(bus, p.parametres);
        } catch (err) {
          alert(err.message);
        }
      });
      d.append(b, el("p", { class: "doux", style: "font-size:12px" },
        "Definitions classiques de Kundur (chap. 4). Le jeu par defaut est celui de son exemple 3.1."));
    })
    .catch((err) => d.append(el("div", { class: "pb erreur" }, err.message)));
  return d;
}

// --- validation -------------------------------------------------------
let minuteurValidation = null;
function validerBientot() {
  clearTimeout(minuteurValidation);
  minuteurValidation = setTimeout(validerMaintenant, 400);
}

async function validerMaintenant() {
  if (!E.tables) return;
  const v = E.version;
  try {
    const d = await api("/api/valider", corpsReseau());
    if (v !== E.version) return;
    const avant = JSON.stringify((E.disj || {}).hors_tension || []);
    E.disj = (E.tables.OUVERTS || []).length
      ? await api("/api/disjoncteurs", corpsReseau()).catch(() => null)
      : null;
    if (v !== E.version) return;
    if (JSON.stringify((E.disj || {}).hors_tension || []) !== avant) dessiner();
    const n = $("#problemes");
    n.textContent = "";
    if (!d.problemes.length) {
      n.append(el("div", { class: "pb ok" }, "Reseau valide."));
      return;
    }
    for (const p of d.problemes) {
      const ligne = el("div", { class: "pb " + p.gravite }, p.message);
      // une reparation s'applique aux tables EDITEES ; en mode zones, la
      // validation porte sur le reseau aplati, qu'on n'edite pas
      if (p.reparation && !E.reseau) {
        ligne.append(" ", el("button", {
          class: "petit reparer",
          title: "Appliquer cette correction (annulable par Ctrl+Z)",
          onclick: () => reparer(p.reparation),
        }, p.reparation.libelle));
      }
      n.append(ligne);
    }
  } catch (err) {
    $("#problemes").textContent = "";
    $("#problemes").append(el("div", { class: "pb erreur" }, err.message));
  }
}

async function reparer(r) {
  try {
    const d = await api("/api/reparer", corpsReseau({ action: r.action, args: r.args }));
    const nom = E.tables.nom;
    delete d.tables.SIMPLIFICATIONS; // celles de l'onglet restent les seules
    Object.assign(E.tables, d.tables);
    E.tables.nom = nom;
    modifie();
    etat("repare : " + r.libelle + " (Ctrl+Z pour annuler)");
  } catch (err) {
    etat("reparation impossible : " + err.message, true);
  }
}

// =====================================================================
// onglet Tables
// =====================================================================
// Les tables et le schema sont DEUX VUES DE LA MEME CHOSE : E.tables. Une
// cellule modifiee ici redessine le schema ; un equipement pose sur le
// schema ajoute ses lignes ici. Cliquer une ligne la montre sur le schema
// de droite ; ce qui est selectionne dans l'onglet Reseau est surligne ici,
// avec les lignes des autres tables qui touchent le meme bus.
const DEFAUT_TEXTE = { T_Loads: { 4: "P" }, T_VSCs: { 3: "PQ" } };

function busSelectionne() {
  if (!E.sel || E.sel.table === "__instance") return null;
  const r = lignes(E.sel.table)[E.sel.index];
  if (!r) return null;
  if (E.sel.table === "Y_DER" || E.sel.table === "T_BusDC") return +r[0];
  return +r[1];
}

function toucheLeBus(table, r, b) {
  if (b === null) return false;
  const cols = COLS(table);
  return cols.some((c, j) => ["BusID", "From", "To", "BusAC", "BusDC"].includes(c) && +r[j] === b);
}

function rendreTables() {
  const n = $("#vue-tables");
  if (!n || !E.tables) return;
  n.textContent = "";
  const meta = el("div", { class: "carte" });
  meta.append(el("h2", {}, "Bases"));
  for (const [k, u] of [["f0", "Hz"], ["Sb", "VA"]]) {
    const lab = el("label", { class: "doux" }, `${k} [${u}]  `);
    const inp = el("input", { type: "number", step: "any", value: E.tables[k] });
    inp.addEventListener("change", () => {
      E.tables[k] = parseFloat(inp.value);
      modifie();
    });
    lab.append(inp);
    meta.append(lab, el("span", {}, "   "));
  }
  n.append(meta);

  const bSel = busSelectionne();
  for (const [table, cols] of Object.entries(E.defauts.tables)) {
    const rows = lignes(table);
    const c = el(rows.length ? "div" : "details", { class: "carte table-edit", style: "margin-top:10px" });
    const titre = el(rows.length ? "h2" : "summary", {}, `${table} -- ${rows.length} ligne(s) `);
    const ajout = el("button", { class: "petit" }, "+ ligne");
    ajout.addEventListener("click", (ev) => {
      ev.preventDefault();
      ajouterLigne(table);
    });
    titre.append(ajout);
    c.append(titre);
    if (!rows.length) {
      c.append(el("p", { class: "doux" }, "Vide. Colonnes : " + cols.join(", ")));
      n.append(c);
      continue;
    }
    const t = el("table");
    const thead = el("tr");
    cols.forEach((x) => thead.append(el("th", {}, x)));
    thead.append(el("th", {}, ""));
    t.append(el("thead", {}, thead));
    const tb = el("tbody");
    rows.forEach((r, i) => {
      const exact = E.sel && E.sel.table === table && E.sel.index === i;
      const tr = el("tr", { class: exact ? "choisi" : toucheLeBus(table, r, bSel) ? "voisin" : "" });
      r.forEach((v, j) => {
        const td = el("td", { class: "chiffre", contenteditable: "true", spellcheck: "false" },
          typeof v === "number" ? String(v) : String(v));
        td.addEventListener("keydown", (ev) => {
          if (ev.key === "Enter") {
            ev.preventDefault();
            td.blur();
          }
        });
        td.addEventListener("blur", () => {
          const brut = td.textContent.trim();
          if (brut === String(v)) return;
          let nv;
          if (typeof v === "number") {
            nv = parseFloat(brut.replace(",", "."));
            if (Number.isNaN(nv)) {
              td.textContent = String(v);
              etat(`${table} : "${brut}" n'est pas un nombre`, true);
              return;
            }
          } else nv = brut;
          const copie = [...r];
          copie[j] = nv;
          rows[i] = copie;
          modifie();
        });
        tr.append(td);
      });
      const sup = el("button", { class: "petit", title: "supprimer la ligne" }, "x");
      sup.addEventListener("click", (ev) => {
        ev.stopPropagation();
        rows.splice(i, 1);
        if (E.sel && E.sel.table === table) E.sel = null;
        modifie();
      });
      tr.append(el("td", {}, sup));
      // cliquer une ligne la selectionne SANS quitter l'onglet : le schema
      // de droite la montre, et les lignes qui touchent le meme bus s'allument
      tr.addEventListener("click", (ev) => {
        if (ev.target.isContentEditable && document.activeElement === ev.target) return;
        E.sel = { table, index: i };
        rendreTables();
        dessiner();
      });
      tb.append(tr);
    });
    t.append(tb);
    c.append(el("div", { class: "defilant" }, t));
    n.append(c);
  }

  const apercu = $("#canevas-tables");
  if (apercu) {
    dessinerSchema(apercu, {});
    $("#tables-sel").textContent = E.sel && E.sel.table !== "__instance" ? `${E.sel.table} [${E.sel.index}]` : "";
  }
}

/** Une ligne de plus : l'identifiant suivant, le reste recopie de la
 *  derniere ligne -- ou des zeros si la table est vide. */
function ajouterLigne(table) {
  const rows = lignes(table);
  const cols = COLS(table);
  let r;
  if (rows.length) r = [...rows[rows.length - 1]];
  else r = cols.map((_, j) => (DEFAUT_TEXTE[table] && DEFAUT_TEXTE[table][j]) || 0);
  r[0] = table === "Y_DER" || table === "T_BusDC" ? nouveauNumero() : nouvelId(table);
  rows.push(r);
  E.sel = { table, index: rows.length - 1 };
  modifie();
}

// le schema de l'onglet Tables se clique aussi : il selectionne la ligne
function initApercuTables() {
  const svg = $("#canevas-tables");
  if (!svg) return;
  svg.addEventListener("pointerdown", (ev) => {
    const s2 = ev.target.closest("[data-sel-table]");
    const b = ev.target.closest("[data-bus]");
    if (s2) E.sel = { table: s2.dataset.selTable, index: +s2.dataset.selIndex };
    else if (b) {
      const table = b.dataset.dc === "1" ? "T_BusDC" : "Y_DER";
      E.sel = { table, index: lignes(table).findIndex((r) => +r[0] === +b.dataset.bus) };
    } else E.sel = null;
    rendreTables();
    dessiner();
  });
}

// =====================================================================
// load flow
// =====================================================================
async function calculerLoadflow() {
  const v = E.version;
  $("#loadflow-etat").textContent = "calcul...";
  try {
    const d = await api("/api/loadflow", corpsReseau());
    if (v !== E.version) return;
    E.res.lf = d;
    $("#loadflow-etat").textContent =
      `${d.iterations} iterations, mismatch ${d.mismatch.toExponential(1)}, ` +
      `pertes ${fmt(d.pertes_P, 4)} pu`;
    rendreLoadflow();
    rendrePoint();
  } catch (err) {
    $("#loadflow-etat").textContent = "";
    $("#table-loadflow").textContent = "";
    $("#table-loadflow").append(el("div", { class: "pb erreur" }, err.message));
  }
}

function rendreLoadflow() {
  const d = E.res.lf;
  if (!d) return;
  const parBus = Object.fromEntries(d.bus.map((b) => [b.id, b]));
  const couleur = (b) => {
    const x = parBus[b];
    if (!x) return "#8a93a5";
    const t = Math.min(1, Math.max(0, (x.V - 0.9) / 0.2)); // 0,9 -> 1,1 pu
    const r = Math.round(30 + 200 * (1 - t)), g = Math.round(90 + 60 * (1 - Math.abs(t - 0.5) * 2)), bl = Math.round(40 + 190 * t);
    return `rgb(${r},${g},${bl})`;
  };
  // En mode zones, le load flow porte sur le reseau APLATI -- ses numeros de
  // bus sont ceux des instances, pas ceux de la vue ouverte. Le dessiner sur
  // la vue colorerait des bus locaux avec les tensions d'autres bus.
  // transits : P entrant par le bus "de", quelle que soit l'orientation de
  // la branche dans la table
  const parBranche = {};
  for (const x of d.branches || []) {
    parBranche[`${x.de}-${x.vers}`] = x.P_de;
    parBranche[`${x.vers}-${x.de}`] = x.P_vers;
  }
  const dessin = () =>
    dessinerSchema($("#canevas-lf"), {
      couleurBus: couleur,
      etiquette: (b) => (parBus[b] ? fmt(parBus[b].V, 3) : ""),
      transit: $("#voir-transits").checked ? (de, vers) => parBranche[`${de}-${vers}`] ?? null : null,
    });
  if (E.reseau && E.plat) avecSource(E.plat, dessin);
  else dessin();

  const n = $("#table-loadflow");
  n.textContent = "";
  const t = el("table");
  const th = el("tr");
  ["bus", "type", "|V| [pu]", "angle [deg]", "P net [MW]", "Q net [Mvar]"].forEach((x) => th.append(el("th", {}, x)));
  t.append(el("thead", {}, th));
  const tb = el("tbody");
  for (const b of d.bus) {
    const tr = el("tr");
    tr.append(el("td", {}, String(b.id)));
    tr.append(el("td", {}, ["", "slack", "PV", "PQ"][b.type] || ""));
    for (const v of [b.V, b.angle]) tr.append(el("td", { class: "chiffre" }, fmt(v, 4)));
    for (const v of [b.P, b.Q]) tr.append(el("td", { class: "chiffre" }, fmt(net0(v * E.tables.Sb / 1e6), 2)));
    tb.append(tr);
  }
  t.append(tb);
  n.append(t);
}

// --- puissance de court-circuit -----------------------------------------
// Une estimation par noeud (g2egrid.court_circuit) : le schema colore chaque
// noeud qui porte des GFL par son SCR, et etiquette les autres de leur Scc.
const couleurSCR = (scr) =>
  scr === null || scr === undefined ? "#667085" : scr < 2 ? "#7a0f0f" : scr < 3 ? "#d92d20" : scr < 5 ? "#dc8a00" : "#00743f";

async function calculerScc() {
  const v = E.version;
  $("#scc-etat").textContent = "calcul...";
  try {
    const d = await api("/api/court_circuit", corpsReseau({
      gfm: $("#scc-gfm").value,
      gfl: $("#scc-gfl").value,
      charges: $("#scc-charges").checked,
      dynamique: $("#scc-dyn").checked,
    }));
    if (v !== E.version) return;
    E.res.scc = d;
    const faibles = d.noeuds.filter((x) => x.SCR !== null && x.SCR < d.seuils.faible).length;
    $("#scc-etat").textContent =
      (d.SCR_min !== null ? `SCR min ${fmt(d.SCR_min, 2)}, SCRIF min ${fmt(d.SCRIF_min, 2)} ; ` : "pas de GFL ; ") +
      `${faibles} noeud(s) a SCR < ${d.seuils.faible}`;
    rendreScc();
  } catch (err) {
    $("#scc-etat").textContent = "";
    $("#table-scc").textContent = "";
    $("#table-scc").append(el("div", { class: "pb erreur" }, err.message));
  }
}

function rendreScc() {
  const d = E.res.scc;
  if (!d) return;
  const parBus = Object.fromEntries(d.noeuds.map((x) => [x.bus, x]));
  const dessin = () =>
    dessinerSchema($("#canevas-scc"), {
      couleurBus: (b) => (parBus[b] ? couleurSCR(parBus[b].SCR) : "#c9ced6"),
      etiquette: (b) => {
        const x = parBus[b];
        if (!x) return "";
        // pas fmt(.., 0) : il retirerait le zero de "40"
        return x.SCR !== null ? `SCR ${fmt(x.SCR, 2)}` : `${Math.round(x.Scc / 1e6)} MVA`;
      },
    });
  if (E.reseau && E.plat) avecSource(E.plat, dessin);
  else dessin();
  const dyn = d.hypotheses.dynamique;
  const n = $("#table-scc");
  n.textContent = "";
  const t = el("table");
  const th = el("tr");
  ["bus", "Un [kV]", "Scc [MVA]", "Ik [kA]", "X/R", "GFL [MVA]", "SCR", "SCRIF", ...(dyn ? ["Scc dq [MVA]"] : [])]
    .forEach((x) => th.append(el("th", {}, x)));
  t.append(el("thead", {}, th));
  const tb = el("tbody");
  for (const x of d.noeuds) {
    const tr = el("tr", { class: x.SCR !== null && x.SCR < d.seuils.faible ? "scc-faible" : null, "data-bus": x.bus });
    const cases = [
      String(x.bus), fmt(x.Un / 1e3, 3), fmt(x.Scc / 1e6, 1), fmt(x.Ik / 1e3, 3), fmt(x.XR, 2),
      x.S_gfl ? fmt(x.S_gfl / 1e6, 2) : "", x.SCR !== null ? fmt(x.SCR, 2) : "",
      x.SCRIF !== null ? fmt(x.SCRIF, 2) : "",
      ...(dyn ? [x.Scc_dyn !== null ? fmt(x.Scc_dyn / 1e6, 1) : "--"] : []),
    ];
    const teinte = { 6: x.SCR, 7: x.SCRIF };
    cases.forEach((c, i) => tr.append(el("td", {
      class: i ? "chiffre" : null,
      style: teinte[i] !== undefined && teinte[i] !== null ? `color:${couleurSCR(teinte[i])};font-weight:600` : null,
    }, c)));
    tb.append(tr);
  }
  t.append(tb);
  n.append(t);
}

// --- le point de fonctionnement ----------------------------------------
// Une ligne par bus qui porte de la production, de la charge, ou un role
// (slack, PV). Tout s'ecrit dans Y_DER : c'est LA table du point de
// fonctionnement, celle que lisent le load flow, les notebooks et le MATLAB.
const TYPES_BUS = { 1: "slack", 2: "PV", 3: "PQ" };

// une puissance de 1e-13 MW est le residu d'arrondi du load flow, pas une
// valeur : l'afficher en notation scientifique fait croire a un transit
const net0 = (v) => (Math.abs(v) < 1e-6 ? 0 : v);

function rendrePoint() {
  const n = $("#table-point");
  if (!n || !E.tables) return;
  n.textContent = "";
  const lignesPoint = lignes("Y_DER")
    .map((r, i) => ({ r, i, eq: equipementsDuBus(+r[0]).filter((e) => !["T_Loads", "T_Shunts"].includes(e.table)) }))
    .filter(({ r, eq }) => eq.length || +r[7] || +r[8] || +r[1] !== 3);
  const lf = E.res.lf ? Object.fromEntries(E.res.lf.bus.map((b) => [b.id, b])) : {};

  const t = el("table");
  const th = el("tr");
  ["bus", "equipement", "type", "V consigne [pu]", "P prod. [MW]", "Q prod. [Mvar]",
   "P charge [MW]", "Q charge [Mvar]", "|V| [pu]", "angle [deg]", "Q prod. calc. [Mvar]"]
    .forEach((x) => th.append(el("th", {}, x)));
  t.append(el("thead", {}, th));
  const tb = el("tbody");
  const Sb = E.tables.Sb;
  for (const { r, i, eq } of lignesPoint) {
    const tr = el("tr");
    tr.append(el("td", {}, String(r[0])));
    tr.append(el("td", { class: "mono" }, eq.map((e) => `${e.famille}${lignes(e.table)[e.index][0]}`).join(", ")));
    const typ = el("select", {});
    for (const [v, nom] of Object.entries(TYPES_BUS)) typ.append(el("option", { value: v, selected: +r[1] === +v ? "" : null }, nom));
    typ.addEventListener("change", () => ecrireY(i, 1, +typ.value));
    tr.append(el("td", {}, typ));
    // chaque consigne n'a de sens que pour certains types : on grise les autres
    const champ = (j, echelle, actif) => {
      const c = el("input", { type: "number", step: "any", value: +(r[j] / echelle).toFixed(6), style: "width:88px" });
      if (!actif) c.setAttribute("disabled", "");
      c.addEventListener("change", () => ecrireY(i, j, parseFloat(c.value) * echelle));
      return el("td", {}, c);
    };
    const ty = +r[1];
    tr.append(champ(3, 1, ty !== 3));
    tr.append(champ(5, 1e6, ty !== 1 && eq.length > 0));
    tr.append(champ(6, 1e6, ty === 3 && eq.length > 0));
    tr.append(champ(7, 1e6, true));
    tr.append(champ(8, 1e6, true));
    const x = lf[+r[0]];
    tr.append(el("td", { class: "chiffre" }, x ? fmt(x.V, 4) : "--"));
    tr.append(el("td", { class: "chiffre" }, x ? fmt(x.angle, 3) : "--"));
    tr.append(el("td", { class: "chiffre" }, x ? fmt(net0((x.Q + x.Ql) * Sb / 1e6), 2) : "--"));
    tb.append(tr);
  }
  t.append(tb);
  n.append(t);

  // le choix du slack : parmi les bus qui portent une unite
  const sel = $("#slack-choisi");
  const garde = sel.value;
  sel.textContent = "";
  for (const { r, eq } of lignesPoint) {
    if (!eq.length) continue;
    sel.append(el("option", { value: r[0], selected: +r[1] === 1 ? "" : null },
      `${eq.map((e) => `${e.famille}${lignes(e.table)[e.index][0]}`).join(", ")} (bus ${r[0]})`));
  }
  if (garde && [...sel.options].some((o) => o.value === garde)) sel.value = garde;
  $("#btn-un-slack").disabled = !!E.reseau;
  $("#btn-un-slack").title = E.reseau
    ? "reseau a zones : regler le point dans la vue du transport, bus par bus"
    : "";
}

function ecrireY(i, j, v) {
  if (Number.isNaN(v)) return;
  const r = [...lignes("Y_DER")[i]];
  r[j] = v;
  lignes("Y_DER")[i] = r;
  modifie();
}

async function unSlackToutPQ() {
  const b = +$("#slack-choisi").value;
  try {
    const d = await api("/api/point/un_slack", corpsReseau({ slack_bus: b }));
    E.tables.Y_DER = d.tables.Y_DER;
    modifie();
    etat(`bus ${b} en slack, tout le reste en PQ : le load flow retombe sur le meme point`);
    await calculerLoadflow();
  } catch (err) {
    etat("impossible : " + err.message, true);
  }
}

// =====================================================================
// plan lineaire commun aux cartes des modes (Modes, Variations, Robustesse)
// =====================================================================
// Toujours une echelle LINEAIRE, des axes epais, des bornes qu'on saisit ou
// qu'on tire au rectangle. Tout est trace ; c'est l'utilisateur qui zoome.
const RENDU_BORNES = {};

function barreBornes(p, redessiner) {
  const n = $("#" + p + "-bornes");
  if (!n) return;
  RENDU_BORNES[p] = redessiner;
  if (n.childElementCount) return;
  const champ = (cle, lab) => {
    const i = el("input", { type: "number", step: "any", id: `${p}-${cle}`, placeholder: "auto", style: "width:92px" });
    i.addEventListener("change", () => RENDU_BORNES[p]());
    return el("div", {}, el("label", { for: `${p}-${cle}` }, lab), i);
  };
  n.append(champ("xmin", "X min"), champ("xmax", "X max"), champ("ymin", "Y min"), champ("ymax", "Y max"));
  const b = el("button", { type: "button", title: "Revenir a l'echelle automatique (ou double-clic sur la carte)" }, "Tout voir");
  b.addEventListener("click", () => {
    ecrireBornes(p, null);
    RENDU_BORNES[p]();
  });
  n.append(b);
}

const CLES_BORNES = [["xmin", "x0"], ["xmax", "x1"], ["ymin", "y0"], ["ymax", "y1"]];

function ecrireBornes(p, b) {
  for (const [k, c] of CLES_BORNES) {
    const i = $(`#${p}-${k}`);
    if (i) i.value = b ? String(+b[c].toPrecision(5)) : "";
  }
}

function minMax(v) {
  let a = Infinity, b = -Infinity;
  for (const x of v) {
    if (x < a) a = x;
    if (x > b) b = x;
  }
  return [a, b];
}

/** Les bornes du trace : celles saisies, l'automatique pour les cases vides. */
function bornesPlan(p, xs, ys, { zeroX = false, zeroY = false } = {}) {
  const auto = (v, zero) => {
    let [a, b] = v.length ? minMax(v) : [-1, 1];
    if (zero) {
      a = Math.min(a, 0);
      b = Math.max(b, 0);
    }
    const pad = 0.05 * (b - a || Math.abs(a) || 1);
    return [a - pad, b + pad];
  };
  const lu = (k) => {
    const i = $(`#${p}-${k}`);
    return i ? parseFloat(i.value) : NaN;
  };
  const [ax0, ax1] = auto(xs, zeroX), [ay0, ay1] = auto(ys, zeroY);
  let x0 = isFinite(lu("xmin")) ? lu("xmin") : ax0, x1 = isFinite(lu("xmax")) ? lu("xmax") : ax1;
  let y0 = isFinite(lu("ymin")) ? lu("ymin") : ay0, y1 = isFinite(lu("ymax")) ? lu("ymax") : ay1;
  if (!(x1 > x0)) [x0, x1] = [ax0, ax1];
  if (!(y1 > y0)) [y0, y1] = [ay0, ay1];
  return { x0, x1, y0, y1 };
}

/** Une graduation en chiffres ordinaires : -300000 plutot que -3.00e+5. */
const graduation = (v) => String(+v.toPrecision(6));

function graduations(a, b, n = 6) {
  const brut = (b - a) / n;
  const e = 10 ** Math.floor(Math.log10(brut));
  const pas = [1, 2, 5, 10].map((k) => k * e).find((k) => brut <= k) || e * 10;
  const out = [];
  for (let v = Math.ceil(a / pas) * pas; v <= b + pas * 1e-9; v += pas) out.push(Math.abs(v) < pas * 1e-9 ? 0 : v);
  return out;
}

/** Grille, axes epais, droites x = 0 et y = 0 ; rend les echelles et un
 *  groupe ecrete a la zone de trace, ou poser les donnees. */
function planLineaire(svg, b, o) {
  const { W, H, ML, MR, MT, MB } = o;
  const PX = (x) => ML + ((x - b.x0) / (b.x1 - b.x0)) * (W - ML - MR);
  const PY = (y) => H - MB - ((y - b.y0) / (b.y1 - b.y0)) * (H - MT - MB);
  const vx = (X) => b.x0 + ((X - ML) / (W - ML - MR)) * (b.x1 - b.x0);
  const vy = (Y) => b.y0 + ((H - MB - Y) / (H - MT - MB)) * (b.y1 - b.y0);
  for (const x of graduations(b.x0, b.x1)) {
    svg.append(el("line", { x1: PX(x), y1: MT, x2: PX(x), y2: H - MB, class: "grille" }));
    svg.append(el("text", { x: PX(x), y: H - MB + 16, "text-anchor": "middle", class: "etiquette" }, graduation(x)));
  }
  for (const y of graduations(b.y0, b.y1, 5)) {
    svg.append(el("line", { x1: ML, y1: PY(y), x2: W - MR, y2: PY(y), class: "grille" }));
    svg.append(el("text", { x: ML - 8, y: PY(y) + 4, "text-anchor": "end", class: "etiquette" }, graduation(y)));
  }
  // les droites x = 0 et y = 0 (l'axe imaginaire est la frontiere de stabilite)
  if (b.x0 < 0 && b.x1 > 0) svg.append(el("line", { x1: PX(0), y1: MT, x2: PX(0), y2: H - MB, class: "axe-zero" }));
  if (b.y0 < 0 && b.y1 > 0) svg.append(el("line", { x1: ML, y1: PY(0), x2: W - MR, y2: PY(0), class: "axe-zero" }));
  // l'abscisse et l'ordonnee, en traits epais
  svg.append(el("line", { x1: ML, y1: H - MB, x2: W - MR, y2: H - MB, class: "axe-fort" }));
  svg.append(el("line", { x1: ML, y1: MT, x2: ML, y2: H - MB, class: "axe-fort" }));
  svg.append(el("text", { x: (ML + W - MR) / 2, y: H - 8, "text-anchor": "middle", class: "etiquette titre-axe" }, o.titreX));
  svg.append(el("text", { x: ML, y: MT - 8, class: "etiquette titre-axe" }, o.titreY));
  const id = svg.id + "-cadre";
  svg.append(el("clipPath", { id }, el("rect", { x: ML, y: MT, width: W - ML - MR, height: H - MT - MB })));
  const g = el("g", { "clip-path": `url(#${id})` });
  svg.append(g);
  return { PX, PY, vx, vy, g, o, b };
}

/** Zoom au rectangle (les bornes se remplissent), double-clic : tout voir.
 *  Les ecouteurs sont poses UNE fois par carte ; le plan courant est relu. */
function zoomPlan(svg, p, P) {
  svg._plan = P;
  svg._bornes = p;
  if (svg._zoomPose) return;
  svg._zoomPose = true;
  let depart = null, rect = null;
  const dans = (q) => {
    const o = svg._plan.o;
    return q[0] >= o.ML && q[0] <= o.W - o.MR && q[1] >= o.MT && q[1] <= o.H - o.MB;
  };
  svg.addEventListener("mousedown", (ev) => {
    const q = pointSVG(svg, ev);
    if (!dans(q)) return;
    depart = q;
    rect = el("rect", { x: q[0], y: q[1], width: 0, height: 0, fill: "#1f4e79", "fill-opacity": 0.12, stroke: "#1f4e79", "stroke-dasharray": "4 3" });
    svg.append(rect);
    ev.preventDefault();
  });
  svg.addEventListener("mousemove", (ev) => {
    if (!depart) return;
    const q = pointSVG(svg, ev);
    rect.setAttribute("x", Math.min(depart[0], q[0]));
    rect.setAttribute("y", Math.min(depart[1], q[1]));
    rect.setAttribute("width", Math.abs(q[0] - depart[0]));
    rect.setAttribute("height", Math.abs(q[1] - depart[1]));
  });
  svg.addEventListener("mouseup", (ev) => {
    if (!depart) return;
    const q = pointSVG(svg, ev), a = depart;
    depart = null;
    if (rect) rect.remove();
    rect = null;
    // un simple clic n'est pas un zoom : il sert a choisir un mode
    if (Math.abs(q[0] - a[0]) < 8 || Math.abs(q[1] - a[1]) < 8) return;
    const P2 = svg._plan;
    ecrireBornes(svg._bornes, {
      x0: P2.vx(Math.min(a[0], q[0])), x1: P2.vx(Math.max(a[0], q[0])),
      y0: P2.vy(Math.max(a[1], q[1])), y1: P2.vy(Math.min(a[1], q[1])),
    });
    RENDU_BORNES[svg._bornes]();
  });
  svg.addEventListener("dblclick", () => {
    ecrireBornes(svg._bornes, null);
    RENDU_BORNES[svg._bornes]();
  });
}

/** Les rayons d'iso-amortissement dans le plan (Re, Im) ou (Re, f). */
function isoAmortissement(P, enFrequence, niveaux = [0.05, 0.1, 0.2]) {
  const { x0, y0, y1 } = P.b;
  for (const zz of niveaux) {
    const k = zz / Math.sqrt(1 - zz * zz); // |Re| = k . |w|
    for (const signe of enFrequence ? [1] : [1, -1]) {
      const yb = signe > 0 ? y1 : y0;
      if (yb * signe <= 0) continue;
      const w = enFrequence ? 2 * Math.PI * yb : yb;
      let xe = -k * Math.abs(w), ye = yb;
      if (xe < x0) { // le rayon sort par la gauche
        xe = x0;
        const wx = Math.abs(x0) / k;
        ye = signe * (enFrequence ? wx / (2 * Math.PI) : wx);
      }
      P.g.append(el("line", { x1: P.PX(0), y1: P.PY(0), x2: P.PX(xe), y2: P.PY(ye), stroke: "#aab2bf", "stroke-dasharray": "4 4" }));
      if (signe > 0) P.g.append(el("text", { x: P.PX(xe) + 4, y: P.PY(ye) + 12, class: "etiquette" }, `${zz * 100} %`));
    }
  }
}

// =====================================================================
// modal
// =====================================================================

// Les grandeurs qu'on peut porter sur les axes de la carte, toujours en
// echelle LINEAIRE : les modes lents se voient en zoomant.
const AXES = {
  re: { nom: "Re(lambda) [1/s]", v: (m) => m.re },
  im: { nom: "Im(lambda) [rad/s]", v: (m) => m.im },
  f: { nom: "frequence amortie [Hz]", v: (m) => m.f },
  f_und: { nom: "frequence propre [Hz]", v: (m) => m.f_und },
  zeta: { nom: "amortissement zeta [%]", v: (m) => m.zeta },
};

async function calculerModal() {
  const v = E.version;
  $("#modal-etat").textContent = "calcul...";
  try {
    const d = await api("/api/t_eig", corpsReseau());
    if (v !== E.version) return;
    E.res.modal = d;
    rendreModal();
  } catch (err) {
    $("#modal-etat").textContent = "";
    $("#table-modes").textContent = "";
    $("#table-modes").append(el("div", { class: "pb erreur" }, err.message));
  }
}

function modesVisibles(d) {
  const osc = $("#osc-seuls").checked;
  const dom = d.domaines ? $("#filtre-domaine").value : "";
  return d.modes.filter((m) => (!osc || m.im > 1e-6) && domaineGarde(m, dom));
}

/** Le filtre par domaine : AC, DC, AC/DC, ou tout mode couple (AC/DC, ou
 *  alternatif porte par deux iles que le continu relie). */
function domaineGarde(m, dom) {
  if (!dom) return true;
  if (dom === "couples") return m.domaine === "AC/DC" || !!m.inter_iles;
  return m.domaine === dom;
}

/** Le marqueur d'un mode : rond (AC), carre (DC), losange (AC/DC). */
function marqueurMode(m, x, y, attrs) {
  if (m.domaine === "DC") return el("rect", { x: x - 5, y: y - 5, width: 10, height: 10, ...attrs });
  if (m.domaine === "AC/DC") return el("polygon", { points: `${x},${y - 7} ${x + 7},${y} ${x},${y + 7} ${x - 7},${y}`, ...attrs });
  return el("circle", { cx: x, cy: y, r: 5.5, ...attrs });
}

function texteDomaine(m) {
  if (!m.domaine) return "";
  const iles = (m.iles || []).map((i) => `${i.nom} ${i.part.toFixed(0)} %`).join(", ");
  return `domaine ${m.domaine} (continu ${m.part_dc.toFixed(0)} %)${m.inter_iles ? ", inter-iles" : ""}` + (iles ? "\n" + iles : "");
}

function rendreModal() {
  const d = E.res.modal;
  if (!d) return;
  const modes = modesVisibles(d);
  const ax = $("#axe-x").value || "re", ay = $("#axe-y").value || "f";
  $("#modal-etat").textContent =
    `${d.n_etats} etats, ${modes.length} modes affiches sur ${d.modes.length}, max Re = ${d.max_re.toExponential(2)} ` +
    (d.stable ? "(stable)" : "(INSTABLE)") +
    (d.jauge ? ` -- ${d.jauge} mode${d.jauge > 1 ? "s" : ""} de jauge ecarte${d.jauge > 1 ? "s" : ""}` : "") +
    (d.domaines ? ` -- AC ${d.domaines.AC}, DC ${d.domaines.DC}, AC/DC ${d.domaines["AC/DC"]}, inter-iles ${d.domaines["inter-iles"]} (* dans la table)` : "");
  $("#bloc-filtre-domaine").hidden = !d.domaines;
  $("#rangee-impact-dc").hidden = !d.domaines;
  $("#modal-etat").style.color = d.stable ? "var(--doux)" : "var(--erreur)";

  const svg = $("#carte-modes");
  svg.textContent = "";
  barreBornes("mod", rendreModal);
  const xs = modes.map((m) => AXES[ax].v(m)), ys = modes.map((m) => AXES[ay].v(m));
  const b = bornesPlan("mod", xs, ys, { zeroX: ax === "re", zeroY: ay === "im" || ay === "re" });
  const P = planLineaire(svg, b, { W: 1000, H: 520, ML: 70, MR: 18, MT: 26, MB: 44, titreX: AXES[ax].nom, titreY: AXES[ay].nom });

  // reperes d'amortissement : les rayons iso-zeta dans le plan (Re, Im) ou
  // (Re, f), le seuil de 5 % des qu'un axe EST l'amortissement
  if (ax === "re" && (ay === "f" || ay === "im")) isoAmortissement(P, ay === "f");
  if (ay === "zeta" && b.y0 < 5 && b.y1 > 5)
    P.g.append(el("line", { x1: P.PX(b.x0), y1: P.PY(5), x2: P.PX(b.x1), y2: P.PY(5), stroke: "#e08b2a", "stroke-dasharray": "5 3" }));
  if (ax === "zeta" && b.x0 < 5 && b.x1 > 5)
    P.g.append(el("line", { x1: P.PX(5), y1: P.PY(b.y0), x2: P.PX(5), y2: P.PY(b.y1), stroke: "#e08b2a", "stroke-dasharray": "5 3" }));

  for (const m of modes) {
    const X = AXES[ax].v(m), Y = AXES[ay].v(m);
    if (X < b.x0 || X > b.x1 || Y < b.y0 || Y > b.y1) continue;
    const c = marqueurMode(m, P.PX(X), P.PY(Y), {
      fill: m.re > 0 ? "#bf2617" : m.zeta < 5 ? "#e08b2a" : "#1f4e79",
      "fill-opacity": 0.85, stroke: m.indice === E.mode ? "#111" : "#fff",
      "stroke-width": m.indice === E.mode ? 2.5 : 1,
    });
    infobulle(c, () =>
      `${fmt(m.re, 3)} ${m.im >= 0 ? "+" : "-"} ${fmt(Math.abs(m.im), 3)}j` + "\n" +
      `f = ${fmt(m.f, 3)} Hz, zeta = ${fmt(m.zeta, 2)} %` + "\n" +
      (m.domaine ? texteDomaine(m) + "\n" : "") +
      m.etats.map((e) => `${e.nom} ${e.part.toFixed(0)} %`).join("\n") + "\n" + "(cliquer pour la forme modale)");
    c.style.cursor = "pointer";
    c.addEventListener("click", () => ouvrirMode(m.indice));
    P.g.append(c);
  }
  $("#legende-modes").textContent =
    (ax === "re" && (ay === "f" || ay === "im") ? "Rayons pointilles : iso-amortissements a 5, 10 et 20 %. " : "") +
    (ax === "zeta" || ay === "zeta" ? "Trait orange : le seuil de 5 % d'amortissement. " : "") +
    "Orange : amorti a moins de 5 % ; rouge : instable." +
    (d.domaines ? " Rond : mode alternatif ; carre : continu ; losange : couple AC/DC." : "");
  zoomPlan(svg, "mod", P);
  svg.dataset.nom = "modes";
  svg.dataset.csv = csvTeig(modes);

  // --- T_eig, le tableau des notebooks --------------------------------
  const n = $("#table-modes");
  n.textContent = "";
  const t = el("table");
  const th = el("tr");
  ["Mode", "Re", "Im", "f_und [Hz]", "f_dam [Hz]", "zeta [%]", ...(d.domaines ? ["Dom.", "DC %"] : []), "Etat 1", "%1", "Etat 2", "%2", "Etat 3", "%3"]
    .forEach((x) => th.append(el("th", {}, x)));
  t.append(el("thead", {}, th));
  const tb = el("tbody");
  for (const m of modes) {
    const tr = el("tr");
    tr.append(el("td", {}, String(m.indice)));
    for (const v of [m.re, m.im, m.f_und, m.f]) tr.append(el("td", { class: "chiffre" }, fmt(v, 4)));
    tr.append(el("td", { class: "chiffre" + (m.zeta < 0 ? " instable" : m.zeta < 5 ? " faible" : "") }, fmt(m.zeta, 2)));
    if (d.domaines) {
      tr.append(el("td", { title: texteDomaine(m) }, (m.domaine || "") + (m.inter_iles ? "*" : "")));
      tr.append(el("td", { class: "chiffre" }, m.part_dc !== undefined ? m.part_dc.toFixed(0) : ""));
    }
    for (let k = 0; k < 3; k++) {
      const e = m.etats[k];
      tr.append(el("td", { class: "mono" }, e ? e.nom : ""));
      tr.append(el("td", { class: "chiffre" }, e ? e.part.toFixed(0) : ""));
    }
    tr.style.cursor = "pointer";
    if (m.indice === E.mode) tr.classList.add("choisi");
    tr.addEventListener("click", () => ouvrirMode(m.indice));
    tb.append(tr);
  }
  t.append(tb);
  n.append(t);
}

async function calculerImpactContinu() {
  const v = E.version;
  $("#impact-dc-etat").textContent = "calcul...";
  const n = $("#impact-resultat");
  n.textContent = "";
  try {
    const d = await api("/api/impact_continu", corpsReseau());
    if (v !== E.version) return;
    E.res.impact_dc = d;
    $("#impact-dc-etat").textContent =
      `${d.n_etats} -> ${d.n_etats_rigide} etats ; max Re ${fmt(d.max_re, 4)} -> ${fmt(d.max_re_rigide, 4)} (continu rigide)`;
    const table = (titre, lignes, avecRigide) => {
      if (!lignes.length) return;
      n.append(el("h3", {}, titre));
      const t = el("table");
      const th = el("tr");
      ["f [Hz]", "zeta [%]", ...(avecRigide ? ["f rigide", "zeta rigide", "ecart [%]"] : []), "part DC [%]", "elements"]
        .forEach((x) => th.append(el("th", {}, x)));
      t.append(el("thead", {}, th));
      const tb = el("tbody");
      for (const m of lignes) {
        const tr = el("tr");
        tr.append(el("td", { class: "chiffre" }, fmt(m.f, 4)));
        tr.append(el("td", { class: "chiffre" + (m.zeta < 0 ? " instable" : m.zeta < 5 ? " faible" : "") }, fmt(m.zeta, 3)));
        if (avecRigide) {
          tr.append(el("td", { class: "chiffre" }, fmt(m.f_rigide, 4)));
          tr.append(el("td", { class: "chiffre" }, fmt(m.zeta_rigide, 3)));
          tr.append(el("td", { class: "chiffre" }, fmt(100 * m.ecart, 3)));
        }
        tr.append(el("td", { class: "chiffre" }, fmt(m.part_dc, 3)));
        tr.append(el("td", { class: "mono" }, m.elements.map((e) => `${e.nom} ${e.part.toFixed(0)} %`).join(", ")));
        tb.append(tr);
      }
      t.append(tb);
      n.append(el("div", { class: "defilant" }, t));
    };
    table("Modes que le continu deplace le plus", d.deplaces, true);
    table("Modes du continu (retires par le continu rigide)", d.retires, false);
    table("Modes perdus avec le continu rigide", d.perdus, false);
    if (d.nouveaux) n.append(el("p", { class: "doux" }, `${d.nouveaux} mode(s) nouveau(x) avec le continu rigide.`));
  } catch (err) {
    $("#impact-dc-etat").textContent = "";
    n.append(el("div", { class: "pb erreur" }, err.message));
  }
}

function csvTeig(modes) {
  return "mode,Re,Im,f_und_Hz,f_dam_Hz,zeta_pct,etat1,pct1,etat2,pct2,etat3,pct3" + "\n" +
    modes.map((m) => [m.indice, m.re, m.im, m.f_und, m.f, m.zeta,
      ...[0, 1, 2].flatMap((k) => (m.etats[k] ? [m.etats[k].nom, m.etats[k].part] : ["", ""]))].join(",")).join("\n");
}

function choisirAxes(x, y) {
  $("#axe-x").value = x;
  $("#axe-y").value = y;
  ecrireBornes("mod", null); // les bornes etaient en unites des ANCIENS axes
  rendreModal();
}

// --- le mode sur le schema ----------------------------------------------
// La forme modale parle en ETATS ; l'exploitant pense en bus et en lignes.
// On projette donc le mode sur le schema unifilaire : chaque bus prend la
// participation cumulee des elements qui y sont raccordes, chaque branche la
// sienne, et une fleche porte la phase de l'etat dominant du bus. Animer fait
// tourner toutes les fleches a la meme vitesse (ralentie) : les ecarts de
// phase, eux, restent ceux du mode.
const TOUR_MODE_S = 4;

function teinteParticipation(x) {
  // gris clair (rien) -> rouge profond (tout), sur une echelle en racine
  const t = Math.sqrt(Math.min(1, Math.max(0, x)));
  const r = Math.round(228 - 37 * t), g = Math.round(231 - 193 * t), b = Math.round(236 - 213 * t);
  return `rgb(${r},${g},${b})`;
}

function arreterAnimationMode() {
  if (E._animMode) cancelAnimationFrame(E._animMode);
  E._animMode = null;
}

function rendreModeSchema() {
  const d = E.res.forme;
  const rangee = $("#rangee-mode-schema");
  arreterAnimationMode();
  if (!d || !d.projection) {
    rangee.hidden = true;
    return;
  }
  rangee.hidden = false;
  const pb = d.projection.bus || {};
  const pmax = Math.max(1e-12, ...Object.values(pb).map((x) => x.part),
    ...(d.projection.branches || []).map((x) => x.part));
  const parBranche = {};
  for (const x of d.projection.branches || []) {
    const cle = [x.de, x.vers].sort((a, b) => a - b).join("-");
    parBranche[cle] = (parBranche[cle] || 0) + x.part;
  }
  $("#mode-schema-titre").textContent = `f = ${fmt(d.f, 3)} Hz, zeta = ${fmt(d.zeta, 2)} %`;
  const svg = $("#schema-mode");
  const dessin = () =>
    dessinerSchema(svg, {
      couleurBus: (b) => (pb[b] ? teinteParticipation(pb[b].part / pmax) : "#e4e7ec"),
      couleurBranche: (br) => {
        const x = parBranche[[br.de, br.vers].sort((a, b) => a - b).join("-")] || 0;
        return { couleur: teinteParticipation(x / pmax), largeur: 2 + 5 * Math.sqrt(x / pmax) };
      },
      etiquette: (b) => (pb[b] && pb[b].part / pmax > 0.02 ? `${fmt(100 * pb[b].part, 1)} %` : ""),
    });
  if (E.reseau && E.plat) avecSource(E.plat, dessin);
  else dessin();

  // les phaseurs, au-dessus du schema
  const g = el("g", { class: "phaseurs" });
  svg.append(g);
  const fleches = [];
  const poser_ = () => {
    for (const [b, x] of Object.entries(pb)) {
      if (x.part / pmax < 0.01) continue;
      const p = versPixels(svg, positionBus(b));
      const L = 14 + 46 * Math.sqrt(x.part / pmax);
      const ligne = el("line", { x1: p[0], y1: p[1], stroke: "#1c2330", "stroke-width": 2.2 });
      const tete = el("polygon", { fill: "#1c2330" });
      infobulle(ligne, `bus ${b} : ${x.etat}` + "\n" + `participation ${fmt(100 * x.part, 2)} %, phase ${fmt(x.phase, 1)} deg`);
      g.append(ligne, tete);
      fleches.push({ p, L, phase: (x.phase * Math.PI) / 180, ligne, tete });
    }
  };
  if (E.reseau && E.plat) avecSource(E.plat, poser_);
  else poser_();
  const placer = (rot) => {
    for (const f of fleches) {
      const a = f.phase + rot;
      const ux = Math.cos(a), uy = -Math.sin(a);
      const x = f.p[0] + f.L * ux, y = f.p[1] + f.L * uy;
      f.ligne.setAttribute("x2", x);
      f.ligne.setAttribute("y2", y);
      f.tete.setAttribute("points",
        `${x + 7 * ux},${y + 7 * uy} ${x - 4 * uy},${y + 4 * ux} ${x + 4 * uy},${y - 4 * ux}`);
    }
  };
  placer(0);
  if (!$("#mode-animer").checked || !fleches.length) return;
  const t0 = performance.now();
  const pas = (t) => {
    if (rangee.hidden) return arreterAnimationMode();
    placer(((t - t0) / 1000) * ((2 * Math.PI) / TOUR_MODE_S));
    E._animMode = requestAnimationFrame(pas);
  };
  E._animMode = requestAnimationFrame(pas);
}

// --- forme modale -----------------------------------------------------
// La table des modes dit QUI oscille ; elle ne dit pas comment les etats se
// placent les uns par rapport aux autres. C'est pourtant la que se lit la
// nature d'un mode : deux vitesses rotor en phase, tout le groupe bouge
// ensemble ; a 180 degres, les machines s'opposent -- un mode inter-machines.
async function ouvrirMode(indice) {
  const v = E.version;
  E.mode = indice;
  $("#rangee-forme").hidden = false;
  $("#forme-titre").textContent = "calcul...";
  try {
    const d = await api("/api/mode", corpsReseau({ indice, n_etats: 10, n_sensibilites: 10 }));
    if (v !== E.version) return;
    E.res.forme = d;
    rendreModeSchema();
    rendreForme();
    rendreModal();  // pour marquer le mode choisi
  } catch (err) {
    $("#forme-titre").textContent = "";
    $("#table-sensibilites").textContent = "";
    $("#table-parametres").textContent = "";
    $("#table-sensibilites").append(el("div", { class: "pb erreur" }, err.message));
  }
}

// Le suffixe d'un etat -- "SM3.Dwr" -> "Dwr" -- sert a ne comparer que ce qui
// est comparable : une vitesse en pu et une puissance en pu n'ont pas le meme
// ordre de grandeur, et tout tracer ensemble ecrase le diagramme.
function suffixe(nom) {
  const i = nom.lastIndexOf(".");
  return i < 0 ? nom : nom.slice(i + 1);
}

function rendreForme() {
  const d = E.res.forme;
  if (!d) return;
  $("#forme-titre").textContent =
    `f = ${fmt(d.f, 3)} Hz, zeta = ${fmt(d.zeta, 2)} %` +
    ` (${fmt(d.re, 3)} ${d.im >= 0 ? "+" : "-"} ${fmt(Math.abs(d.im), 3)}j)`;

  // le choix des suffixes se refait a chaque mode : ils changent avec lui
  const sel = $("#forme-suffixe");
  const garde = sel.value;
  const sufs = [...new Set(d.forme.map((e) => suffixe(e.nom)))];
  sel.textContent = "";
  sel.append(el("option", { value: "" }, "tous"));
  for (const x of sufs) sel.append(el("option", { value: x }, x));
  sel.value = sufs.includes(garde) ? garde : "";

  const filtre = sel.value;
  const pts = d.forme.filter((e) => !filtre || suffixe(e.nom) === filtre);

  const svg = $("#polaire");
  svg.textContent = "";
  const W = 520, H = 420, cx = W / 2, cy = H / 2 - 6, Rmax = 150;
  const mmax = Math.max(1e-30, ...pts.map((e) => e.module));

  for (const r of [0.25, 0.5, 0.75, 1]) {
    svg.append(el("circle", { cx, cy, r: r * Rmax, fill: "none", class: "grille" }));
    svg.append(el("text", { x: cx + 3, y: cy - r * Rmax + 11, class: "etiquette" }, r.toFixed(2)));
  }
  for (let a = 0; a < 360; a += 30) {
    const t = (a * Math.PI) / 180;
    svg.append(el("line", {
      x1: cx, y1: cy, x2: cx + Rmax * Math.cos(t), y2: cy - Rmax * Math.sin(t),
      class: a % 90 === 0 ? "axe" : "grille",
    }));
  }
  for (const [a, txt] of [[0, "0"], [90, "90"], [180, "180"], [270, "-90"]]) {
    const t = (a * Math.PI) / 180;
    svg.append(el("text", {
      x: cx + (Rmax + 16) * Math.cos(t), y: cy - (Rmax + 16) * Math.sin(t) + 4,
      "text-anchor": "middle", class: "etiquette",
    }, txt + " deg"));
  }

  // une couleur par element : c'est l'element qu'on compare, pas l'etat
  const elems = [...new Set(pts.map((e) => e.nom.split(".")[0]))];
  const teintes = ["#1f4e79", "#bf2617", "#2a7a4b", "#e08b2a", "#6b3fa0", "#0f7d8c"];

  for (const e of pts) {
    const t = (e.phase * Math.PI) / 180;
    const r = (e.module / mmax) * Rmax;
    const x = cx + r * Math.cos(t), y = cy - r * Math.sin(t);
    const col = teintes[elems.indexOf(e.nom.split(".")[0]) % teintes.length];
    svg.append(el("line", { x1: cx, y1: cy, x2: x, y2: y, stroke: col, "stroke-width": 2 }));
    const c = el("circle", { cx: x, cy: y, r: 4.5, fill: col, stroke: "#fff", "stroke-width": 1 });
    infobulle(c, `${e.nom}` + "\n" + `module ${fmt(e.module, 4)} (relatif ${fmt(e.module / mmax, 3)})` + "\n" + `phase ${fmt(e.phase, 1)} deg`);
    svg.append(c);
    if (r > 0.18 * Rmax) {
      svg.append(el("text", {
        x: cx + (r + 12) * Math.cos(t), y: cy - (r + 12) * Math.sin(t) + 4,
        "text-anchor": Math.cos(t) < -0.3 ? "end" : Math.cos(t) > 0.3 ? "start" : "middle",
        class: "etiquette", fill: col,
      }, e.nom));
    }
  }

  svg.dataset.nom = "forme_modale";
  svg.dataset.csv =
    "etat,module,module_relatif,phase_deg" + "\n" +
    d.forme.map((e) => [e.nom, e.module, e.module / mmax, e.phase].join(",")).join("\n");

  const n = $("#table-sensibilites");
  n.textContent = "";
  const tb = el("tbody");
  for (const e of d.sensibilites) {
    const tr = el("tr");
    tr.append(el("td", { class: "mono" }, e.ligne));
    tr.append(el("td", { class: "mono" }, e.colonne));
    tr.append(el("td", { class: "chiffre" }, fmt(e.valeur, 4)));
    tb.append(tr);
  }
  const t = el("table");
  const th = el("tr");
  ["d/dt de...", "...par rapport a", "|d lambda / d A|"].forEach((x) => th.append(el("th", {}, x)));
  t.append(el("thead", {}, th), tb);
  n.append(t);

  // le PARAMETRE qui deplace le mode, et plus seulement le coefficient de A
  const np_ = $("#table-parametres");
  np_.textContent = "";
  if (d.parametres_erreur) {
    np_.append(el("div", { class: "pb" }, d.parametres_erreur));
    return;
  }
  const tbp = el("tbody");
  for (const e of d.parametres || []) {
    const tr = el("tr");
    tr.append(el("td", { class: "mono" }, `${e.element}.${e.parametre}`));
    tr.append(el("td", { class: "chiffre" }, fmt(e.valeur, 4)));
    tr.append(el("td", { class: "chiffre" }, fmt(e.re, 4)));
    tr.append(el("td", { class: "chiffre" }, fmt(e.df, 4)));
    tr.append(el("td", { class: "chiffre" }, fmt(e.dzeta, 3)));
    tbp.append(tr);
  }
  const tp = el("table");
  const thp = el("tr");
  ["parametre", "valeur", "p.d Re", "p.d f [Hz]", "p.d zeta [%]"].forEach((x) => thp.append(el("th", {}, x)));
  tp.append(el("thead", {}, thp), tbp);
  np_.append(tp);
}

// =====================================================================
// temporel
// =====================================================================
async function chargerSignaux() {
  const v = E.version;
  const d = await api("/api/signaux", corpsReseau());
  if (v !== E.version) return;
  E.res.signaux = d;
  const e = $("#entree");
  const garde = e.value;
  e.textContent = "";
  for (const s of d.entrees) e.append(el("option", { value: s }, s));
  // w_ref vient en tete de la liste parce que c'est l'entree GLOBALE, la
  // pulsation du repere -- et c'est le pire premier echelon qu'on puisse
  // proposer : 0,01 pu dessus fait sauter le repere de tout le reseau d'un
  // demi-hertz, ce qui ecarte les deux voies de 386 % sur smib et laisse
  // croire que le modele lineaire est faux. On preselectionne une vraie
  // consigne d'equipement, celle qu'on perturbe pour de bon.
  const consigne =
    d.entrees.find((s) => s.endsWith(".p_ref")) ||
    d.entrees.find((s) => s.endsWith(".v_ref")) ||
    d.entrees.find((s) => s.includes(".")) ||
    d.entrees[0];
  e.value = d.entrees.includes(garde) && garde ? garde : consigne;
  const x = $("#etats");
  x.textContent = "";
  for (const s of d.etats) x.append(el("option", { value: s }, s));
  // une preselection utile : les vitesses de rotor, sinon les premiers etats
  const pref = d.etats.filter((s) => s.endsWith(".Dwr") || s.endsWith(".w") || s.endsWith(".pm")).slice(0, 3);
  for (const o of x.options) o.selected = pref.length ? pref.includes(o.value) : o.index < 2;

  // Les mesures sont ce qu'on compare a un enregistrement : un module de
  // tension, une frequence, un transit. On en preselectionne une paire
  // parlante -- la frequence filtree d'un bus et sa tension -- plutot que
  // rien, pour que le premier clic sur Simuler montre deja quelque chose.
  const q = $("#mesures");
  const gardeQ = Array.from(q.selectedOptions).map((o) => o.value);
  q.textContent = "";
  for (const s2 of d.mesures || []) q.append(el("option", { value: s2 }, s2));
  const prefQ = gardeQ.length
    ? gardeQ
    : [(d.mesures || []).find((s2) => s2.startsWith("fm_")), (d.mesures || []).find((s2) => s2.startsWith("V_"))].filter(Boolean);
  for (const o of q.options) o.selected = prefQ.includes(o.value);

  rendreScenario();
}

// --- scenario ---------------------------------------------------------
// Un scenario est une liste d'evenements dates. Le cas courant reste
// l'echelon unique du formulaire ; cette table ne le remplace que si on
// coche la case, pour ne pas alourdir le geste le plus frequent.
// Deux sortes d'evenements dans le meme scenario : une CONSIGNE qui saute ou
// monte en rampe (echelon, rampe), et une MANOEUVRE qui change le reseau
// (ouverture d'un disjoncteur, variation d'une charge, saut de phase d'un bus).
const GENRES_CONSIGNE = ["echelon", "rampe"];
const GENRES_MANOEUVRE = ["ouverture", "charge", "saut_de_phase", "defaut", "elimination", "defaut_dc", "elimination_dc", "blocage"];

function ciblesManoeuvre(genre) {
  const m = (E.res.signaux && E.res.signaux.manoeuvres) || { elements: [], bus: [], bus_dc: [] };
  if (["saut_de_phase", "defaut", "elimination"].includes(genre)) return m.bus.map((b) => "bus" + b);
  if (["defaut_dc", "elimination_dc"].includes(genre)) return (m.bus_dc || []).map((b) => "bus" + b);
  if (genre === "blocage") return m.elements.filter((e) => e.genre === "convertisseur").map((e) => e.nom);
  if (genre === "charge") return m.elements.filter((e) => e.genre === "charge").map((e) => e.nom);
  return m.elements.map((e) => e.nom);
}

function rendreScenario() {
  const tb = $("#table-scenario tbody");
  if (!tb) return;
  rendreFrise();
  tb.textContent = "";
  const entrees = (E.res.signaux && E.res.signaux.entrees) || [];
  E.scenario.forEach((ev, i) => {
    const tr = el("tr");
    const champ = (cle, pas, largeur = 80) => {
      const n = el("input", { type: "number", value: ev[cle] ?? 0, step: pas, style: `width:${largeur}px` });
      n.addEventListener("change", () => {
        ev[cle] = parseFloat(n.value);
        if (cle === "instant") rendreFrise();
      });
      return n;
    };
    tr.append(el("td", {}, champ("instant", 0.05)));
    const g = el("select");
    for (const x of [...GENRES_CONSIGNE, ...GENRES_MANOEUVRE]) g.append(el("option", { value: x }, x.replaceAll("_", " ")));
    g.value = ev.genre;
    g.addEventListener("change", () => {
      ev.genre = g.value;
      if (GENRES_MANOEUVRE.includes(ev.genre)) ev.cible = ciblesManoeuvre(ev.genre)[0] || "";
      // impedance de defaut par defaut : alternative, ou continue pole a pole
      if (ev.genre === "defaut") Object.assign(ev, { R: 1e-3, X: 1e-2 });
      if (ev.genre === "defaut_dc") Object.assign(ev, { R: 1e-4, X: 1e-3 });
      rendreScenario();
    });
    tr.append(el("td", {}, g));
    const sel = el("select", { style: "max-width:190px" });
    const manoeuvre = GENRES_MANOEUVRE.includes(ev.genre);
    for (const x of manoeuvre ? ciblesManoeuvre(ev.genre) : entrees) sel.append(el("option", { value: x }, x));
    sel.value = manoeuvre ? ev.cible : ev.entree;
    sel.addEventListener("change", () => (manoeuvre ? (ev.cible = sel.value) : (ev.entree = sel.value)));
    tr.append(el("td", {}, sel));
    const params = el("td", {});
    if (ev.genre === "echelon") params.append("amplitude ", champ("amplitude", 0.005));
    else if (ev.genre === "rampe") params.append("amplitude ", champ("amplitude", 0.005), " sur ", champ("duree", 0.05, 60), " s");
    else if (ev.genre === "charge") params.append("dP ", champ("dp", 0.05, 60), " dQ ", champ("dq", 0.05, 60), " (relatif)");
    else if (ev.genre === "saut_de_phase") params.append(champ("angle", 1, 60), " deg");
    else if (ev.genre === "defaut") params.append("R ", champ("R", 0.001, 70), " X ", champ("X", 0.001, 70), " pu");
    else if (ev.genre === "elimination") params.append(el("span", { class: "doux" }, "le defaut est retire"));
    else if (ev.genre === "defaut_dc") params.append("R ", champ("R", 0.0001, 70), " X ", champ("X", 0.0001, 70), " pu (pole a pole)");
    else if (ev.genre === "elimination_dc") params.append(el("span", { class: "doux" }, "le defaut continu est retire"));
    else if (ev.genre === "blocage") params.append(el("span", { class: "doux" }, "le convertisseur se bloque"));
    else params.append(el("span", { class: "doux" }, "le disjoncteur s'ouvre"));
    tr.append(params);
    const sup = el("button", {}, "x");
    sup.addEventListener("click", () => { E.scenario.splice(i, 1); rendreScenario(); });
    tr.append(el("td", {}, sup));
    tb.append(tr);
  });
}

// --- frise du scenario --------------------------------------------------
const COULEUR_EVT = {
  echelon: "#1f4e79", rampe: "#0f7d8c", ouverture: "#bf2617", charge: "#e08b2a",
  saut_de_phase: "#6b3fa0", defaut: "#8b0000", elimination: "#2a7a4b",
  defaut_dc: "#5b1a7a", elimination_dc: "#2a7a4b", blocage: "#444444",
};

function bornesFrise() {
  const fin = Math.max(parseFloat($("#duree").value) || 1, ...E.scenario.map((e) => (e.instant || 0) + 0.05));
  return [0, fin];
}

function rendreFrise() {
  const svg = $("#frise");
  if (!svg) return;
  svg.textContent = "";
  const [t0, t1] = bornesFrise();
  const X = (t) => 30 + ((t - t0) / (t1 - t0)) * 940;
  svg.append(el("line", { x1: 30, y1: 60, x2: 970, y2: 60, stroke: "#344054", "stroke-width": 2 }));
  for (const t of graduations(t0, t1, 8)) {
    svg.append(el("line", { x1: X(t), y1: 56, x2: X(t), y2: 64, stroke: "#344054" }));
    svg.append(el("text", { x: X(t), y: 80, "text-anchor": "middle", class: "etiquette" }, graduation(t) + " s"));
  }
  E.scenario.forEach((ev, i) => {
    const x = X(ev.instant || 0);
    const c = COULEUR_EVT[ev.genre] || "#344054";
    const g = el("g", { style: "cursor:ew-resize", "data-i": i });
    g.append(el("line", { x1: x, y1: 22, x2: x, y2: 60, stroke: c, "stroke-width": 2 }));
    g.append(el("polygon", { points: `${x},14 ${x + 8},22 ${x},30 ${x - 8},22`, fill: c, stroke: "#fff" }));
    const quoi = GENRES_MANOEUVRE.includes(ev.genre) ? ev.cible : ev.entree;
    g.append(el("text", { x: x + 10, y: 18, class: "etiquette", fill: c }, `${ev.genre.replaceAll("_", " ")} ${quoi || ""}`));
    infobulle(g, `${fmt(ev.instant, 4)} s : ${ev.genre} ${quoi || ""}`);
    g.addEventListener("pointerdown", (pe) => {
      pe.stopPropagation();
      svg.setPointerCapture(pe.pointerId);
      E._friseDrag = { i, x0: pe.clientX, t0: ev.instant || 0, bouge: false };
    });
    svg.append(g);
  });
  // le curseur du rejeu
  if (E._rejeuT !== undefined && E._rejeuT !== null) {
    const x = X(E._rejeuT);
    svg.append(el("line", { x1: x, y1: 4, x2: x, y2: 70, stroke: "#1c2330", "stroke-width": 1.5, "stroke-dasharray": "3 2" }));
  }
}

function initFrise() {
  const svg = $("#frise");
  const versT = (clientX) => {
    const r = svg.getBoundingClientRect();
    const [t0, t1] = bornesFrise();
    const u = ((clientX - r.left) / r.width) * 1000;
    return Math.max(0, t0 + ((u - 30) / 940) * (t1 - t0));
  };
  svg.addEventListener("pointermove", (pe) => {
    const d = E._friseDrag;
    if (!d) return;
    const r = svg.getBoundingClientRect();
    const [t0, t1] = bornesFrise();
    const dt = ((pe.clientX - d.x0) / r.width) * 1000 / 940 * (t1 - t0);
    if (Math.abs(pe.clientX - d.x0) > 2) d.bouge = true;
    E.scenario[d.i].instant = Math.max(0, +(d.t0 + dt).toFixed(4));
    rendreFrise();
  });
  svg.addEventListener("pointerup", () => {
    const d = E._friseDrag;
    E._friseDrag = null;
    if (!d) return;
    if (!d.bouge) {
      const ligne = $("#table-scenario tbody").children[d.i];
      if (ligne) {
        ligne.scrollIntoView({ block: "nearest" });
        ligne.classList.add("eclair");
        setTimeout(() => ligne.classList.remove("eclair"), 900);
      }
      return;
    }
    E.scenario.sort((a, b) => (a.instant || 0) - (b.instant || 0));
    rendreScenario();
  });
  svg.addEventListener("dblclick", (pe) => {
    ajouterEvenement();
    E.scenario[E.scenario.length - 1].instant = +versT(pe.clientX).toFixed(3);
    E.scenario.sort((a, b) => (a.instant || 0) - (b.instant || 0));
    rendreScenario();
  });
}

// --- rejeu sur le schema --------------------------------------------------
function couleurTension(V) {
  if (V === null || V === undefined || !Number.isFinite(V)) return "#c9ced6";
  const t = Math.min(1, Math.max(0, (V - 0.9) / 0.2)); // 0,9 -> 1,1 pu
  const r = Math.round(30 + 200 * (1 - t)), g = Math.round(90 + 60 * (1 - Math.abs(t - 0.5) * 2)), b = Math.round(40 + 190 * t);
  return `rgb(${r},${g},${b})`;
}

function preparerRejeu() {
  const d = E.res.temporel;
  const V = Object.keys((d && d.mesures) || {}).filter((k) => k.startsWith("V_bus"));
  $("#rangee-rejeu").hidden = !V.length;
  if (!V.length) return;
  const curseur = $("#rejeu-t");
  curseur.max = d.t.length - 1;
  curseur.value = 0;
  arreterRejeu();
  rendreRejeu(0);
}

function rendreRejeu(i) {
  const d = E.res.temporel;
  if (!d) return;
  const t = d.t[i];
  E._rejeuT = t;
  $("#rejeu-titre").textContent = `t = ${fmt(t, 4)} s`;
  const V = (b) => {
    const s = d.mesures["V_bus" + b];
    return s ? s[i] : null;
  };
  const dessin = () =>
    dessinerSchema($("#schema-rejeu"), {
      couleurBus: (b) => couleurTension(V(b)),
      etiquette: (b) => (V(b) !== null && V(b) !== undefined ? fmt(V(b), 3) : ""),
    });
  if (E.reseau && E.plat) avecSource(E.plat, dessin);
  else dessin();
  const passes = (d.scenario || []).filter((e) => (e.instant || 0) <= t);
  $("#rejeu-evts").textContent = passes.length
    ? "evenements passes : " + passes.map((e) => `${fmt(e.instant, 3)} s ${e.genre.replaceAll("_", " ")} ${GENRES_MANOEUVRE.includes(e.genre) ? e.cible : e.entree}`).join(" ; ")
    : "";
  rendreFrise();
}

function arreterRejeu() {
  if (E._rejeuAnim) cancelAnimationFrame(E._rejeuAnim);
  E._rejeuAnim = null;
  const b = $("#btn-rejeu");
  if (b) b.textContent = "Lecture";
}

function basculerRejeu() {
  if (E._rejeuAnim) return arreterRejeu();
  const d = E.res.temporel;
  if (!d) return;
  const curseur = $("#rejeu-t");
  const n = d.t.length;
  if (+curseur.value >= n - 1) curseur.value = 0;
  const duree = (parseFloat($("#rejeu-duree").value) || 8) * 1000;
  const i0 = +curseur.value, debut = performance.now();
  $("#btn-rejeu").textContent = "Pause";
  const pas = (maintenant) => {
    const i = Math.min(n - 1, i0 + Math.floor(((maintenant - debut) / duree) * (n - 1)));
    curseur.value = i;
    rendreRejeu(i);
    if (i >= n - 1) return arreterRejeu();
    E._rejeuAnim = requestAnimationFrame(pas);
  };
  E._rejeuAnim = requestAnimationFrame(pas);
}

function ajouterEvenement() {
  const entrees = (E.res.signaux && E.res.signaux.entrees) || [];
  const dernier = E.scenario[E.scenario.length - 1];
  E.scenario.push({
    instant: dernier ? dernier.instant + 0.5 : 0.1,
    entree: dernier && dernier.entree ? dernier.entree : $("#entree").value || entrees[0] || "",
    amplitude: dernier && dernier.amplitude ? -dernier.amplitude : parseFloat($("#amplitude").value) || 0.01,
    genre: "echelon",
    duree: 0,
    dp: 0.1,
    dq: 0.1,
    angle: 10,
    R: 0.001,
    X: 0.01,
    cible: "",
  });
  $("#scenario-actif").checked = true;
  $("#bloc-scenario").open = true;
  rendreScenario();
}

async function simuler() {
  const v = E.version;
  $("#temporel-etat").textContent = "simulation...";
  const etats = Array.from($("#etats").selectedOptions).map((o) => o.value);
  const quoi = Array.from($("#mesures").selectedOptions).map((o) => o.value);
  const scen = $("#scenario-actif").checked && E.scenario.length ? E.scenario : null;
  // le rejeu demande la tension de chaque bus ; seules les mesures choisies
  // s'affichent en courbes
  const tousV = $("#rejeu-actif").checked
    ? ((E.res.signaux && E.res.signaux.mesures) || []).filter((m) => m.startsWith("V_"))
    : [];
  const demandees = [...new Set([...quoi, ...tousV])];
  try {
    const d = await api("/api/temporel", corpsReseau({
      entree: $("#entree").value,
      amplitude: parseFloat($("#amplitude").value),
      t0: parseFloat($("#t0").value),
      duree: parseFloat($("#duree").value),
      n: 401,
      etats: etats.length ? etats : null,
      voie: $("#voie").value,
      evenements: scen,
      mesures: demandees.length ? demandees : null,
      methode: $("#integrateur").value,
      pas: parseFloat($("#pas").value) ? parseFloat($("#pas").value) * 1e-6 : null,
      limiteurs: $("#limiteurs").checked,
      protections_dc: !!($("#protections-dc") && $("#protections-dc").checked),
    }));
    if (v !== E.version) return;
    const omise = $("#voie").value !== "non_lineaire" && !d.voies.includes("lineaire");
    $("#temporel-etat").textContent = d.ecart
      ? `ecart lineaire / non lineaire : ${d.ecart.absolu.toExponential(2)} (${(100 * d.ecart.relatif).toFixed(2)} % de l'amplitude)`
      : omise
        ? "voie lineaire omise : elle ne suit qu'un echelon unique"
        : "";
    E.res.temporel = { ...d, affichees: quoi, scenario: scen ? JSON.parse(JSON.stringify(scen)) : [] };
    rendreTemporel({ ...d, mesures: Object.fromEntries(Object.entries(d.mesures || {}).filter(([k]) => quoi.includes(k) || !tousV.includes(k))) });
    preparerRejeu();
  } catch (err) {
    $("#temporel-etat").textContent = "";
    const n = $("#graphes-temporel");
    n.textContent = "";
    n.append(el("div", { class: "pb erreur" }, err.message));
  }
}

function rendreTemporel(d) {
  const n = $("#graphes-temporel");
  n.textContent = "";
  // les limiteurs qui ont joue : quand, qui, dans quel sens
  const b = $("#basculements");
  b.textContent = "";
  if (d.integrateur) b.append(el("span", {}, `integrateur : ${d.integrateur}. `));
  const bs = d.basculements || [];
  if (bs.length) {
    const nom = { mode_lim: "limitation de courant", mode_e2: "butee d'excitation" };
    const lignes_ = bs.slice(0, 12).map((x) =>
      `${fmt(x.t, 4)} s : ${x.element} ${x.valeur ? "entre en" : "sort de"} ${nom[x.drapeau] || x.drapeau}`);
    b.append(el("span", {}, `${bs.length} basculement(s) de limiteur -- ` + lignes_.join(" ; ") + (bs.length > 12 ? " ..." : "")));
  }
  // les protections qui ont declenche, et ce qu'elles ont manoeuvre
  const ds = d.declenchements || [];
  if (ds.length) {
    const lignes_ = ds.map((x) =>
      `${fmt(1e3 * x.t, 4)} ms ${x.protection.replace(":", " ")} -> ` +
      x.manoeuvres.map((m) => `${m.genre} ${m.cible}`).join(", "));
    b.append(el("div", { class: "declenchements" }, `${ds.length} declenchement(s) de protection : ` + lignes_.join(" ; ")));
  }
  // les mesures d'abord : c'est ce qu'on compare a un enregistrement, les
  // etats bruts viennent apres pour comprendre ce qui s'est passe dedans
  for (const [nom, v] of Object.entries(d.mesures || {})) {
    const c = el("div", { class: "carte" });
    c.append(el("h2", {}, nom + unite(nom)));
    const g = graphe(d.t, { mesure: v }, nom);
    g.dataset.nom = nom;
    g.dataset.csv = "t," + nom + "\n" + d.t.map((t, i) => t + "," + v[i]).join("\n");
    c.append(g);
    n.append(c);
  }
  for (const [nom, series] of Object.entries(d.series)) {
    const c = el("div", { class: "carte" });
    c.append(el("h2", {}, nom));
    const g = graphe(d.t, series, nom);
    const voies = Object.keys(series);
    g.dataset.nom = nom;
    g.dataset.csv =
      ["t"].concat(voies).join(",") +
      "\n" +
      d.t.map((t, i) => [t].concat(voies.map((v) => series[v][i])).join(",")).join("\n");
    c.append(g);
    n.append(c);
  }
}

const STYLE_VOIE = {
  lineaire: { stroke: "#1f4e79", "stroke-width": 1.8 },
  non_lineaire: { stroke: "#bf2617", "stroke-width": 1.6, "stroke-dasharray": "5 3" },
  mesure: { stroke: "#2a7a4b", "stroke-width": 1.8 },
};

// Une mesure porte son unite dans son prefixe ; l'afficher evite de lire un
// module de tension comme une frequence.
function unite(nom) {
  if (nom.startsWith("V_")) return " [pu]";
  if (nom.startsWith("angle_")) return " [deg]";
  if (nom.startsWith("f_") || nom.startsWith("fm_")) return " [Hz]";
  if (nom.startsWith("P_") || nom.startsWith("Q_")) return " [pu]";
  return "";
}

function graphe(t, series, titre) {
  const W = 560, H = 260, ML = 64, MR = 14, MT = 12, MB = 32;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, class: "graphe", style: "width:100%" });
  // un NaN rendu null (grandeur d'un element retire par une manoeuvre) est
  // un TROU dans la courbe, pas un zero : on l'ecarte des bornes et on coupe
  // le trace a cet endroit
  const vals = Object.values(series).flat().filter((v) => v !== null && Number.isFinite(v));
  let y0 = Math.min(...vals), y1 = Math.max(...vals);
  if (!(isFinite(y0) && isFinite(y1))) { y0 = -1; y1 = 1; }
  if (y1 - y0 < 1e-12) { y0 -= 1e-9; y1 += 1e-9; }
  const m = (y1 - y0) * 0.08;
  y0 -= m; y1 += m;
  const px = (x) => ML + ((x - t[0]) / (t[t.length - 1] - t[0] || 1)) * (W - ML - MR);
  const py = (y) => H - MB - ((y - y0) / (y1 - y0)) * (H - MT - MB);

  for (let k = 0; k <= 4; k++) {
    const y = y0 + ((y1 - y0) * k) / 4;
    svg.append(el("line", { x1: ML, y1: py(y), x2: W - MR, y2: py(y), class: "grille" }));
    svg.append(el("text", { x: ML - 6, y: py(y) + 4, "text-anchor": "end", class: "etiquette" }, fmt(y, 4)));
  }
  for (let k = 0; k <= 4; k++) {
    const x = t[0] + ((t[t.length - 1] - t[0]) * k) / 4;
    svg.append(el("text", { x: px(x), y: H - MB + 14, "text-anchor": "middle", class: "etiquette" }, fmt(x, 3)));
  }
  if (y0 < 0 && y1 > 0) svg.append(el("line", { x1: ML, y1: py(0), x2: W - MR, y2: py(0), class: "axe" }));
  svg.append(el("text", { x: (ML + W) / 2, y: H - 4, "text-anchor": "middle", class: "etiquette" }, "temps [s]"));

  for (const [voie, y] of Object.entries(series)) {
    let morceau = [];
    const poser = () => {
      if (morceau.length > 1)
        svg.append(el("polyline", Object.assign({ points: morceau.join(" "), fill: "none" }, STYLE_VOIE[voie] || {})));
      morceau = [];
    };
    y.forEach((v, i) => {
      if (v === null || !Number.isFinite(v)) poser();
      else morceau.push(`${px(t[i])},${py(v)}`);
    });
    poser();
  }
  return svg;
}

// =====================================================================
// export des figures : PNG, SVG, CSV
// =====================================================================
/* Un observateur de mutations decore CHAQUE figure au fur et a mesure
 * qu'elle apparait, au lieu que chaque onglet ait a s'y inscrire : ajouter
 * un graphe ailleurs dans la page lui donne sa barre d'export sans une
 * ligne de plus. */
function barreExport(svg) {
  if (svg.dataset.exportable) return;
  svg.dataset.exportable = "1";
  const barre = el("div", { class: "export-figure" });
  const ajoute = (nom, fn) => barre.append(el("button", { onclick: fn }, nom));
  ajoute("SVG", () => telecharger(sourceSVG(svg), "image/svg+xml", "figure.svg"));
  ajoute("PNG", () => versPNG(svg));
  // le CSV se lit AU CLIC : la barre est posee des que la figure apparait,
  // souvent avant que le calcul ait rempli dataset.csv
  ajoute("CSV", () => {
    const csv = svg.dataset.csv;
    if (csv) telecharger(csv, "text/csv", (svg.dataset.nom || "donnees") + ".csv");
  });
  svg.parentNode.insertBefore(barre, svg);
}

/** Le SVG isole de sa page a besoin de ses styles EN LIGNE : les regles du
 *  site sont des selecteurs descendants, qui cessent de s'appliquer des que
 *  la figure quitte le document. */
function sourceSVG(svg) {
  const copie = svg.cloneNode(true);
  const src = svg.querySelectorAll("*");
  copie.querySelectorAll("*").forEach((n, i) => {
    const c = getComputedStyle(src[i]);
    // font-weight et stroke-linecap : sans eux, le "GFM" gras des symboles
    // et les bouts arrondis des traits se perdent a l'export
    for (const prop of ["fill", "stroke", "stroke-width", "stroke-dasharray", "stroke-linecap",
                        "font-size", "font-family", "font-weight", "fill-opacity"])
      n.style[prop] = c.getPropertyValue(prop);
  });
  copie.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  copie.setAttribute("style", "background:#fff");
  return new XMLSerializer().serializeToString(copie);
}

function versPNG(svg) {
  const vb = svg.viewBox.baseVal;
  const img = new Image();
  img.onload = () => {
    const c = document.createElement("canvas");
    c.width = vb.width * 2;
    c.height = vb.height * 2;
    const ctx = c.getContext("2d");
    ctx.fillStyle = "#fff";
    ctx.fillRect(0, 0, c.width, c.height);
    ctx.drawImage(img, 0, 0, c.width, c.height);
    c.toBlob((b) => {
      const a = el("a", { href: URL.createObjectURL(b), download: "figure.png" });
      document.body.append(a);
      a.click();
      a.remove();
    });
  };
  img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(sourceSVG(svg));
}

function telecharger(texte, type, nom) {
  const b = new Blob([texte], { type });
  const a = el("a", { href: URL.createObjectURL(b), download: nom });
  document.body.append(a);
  a.click();
  a.remove();
}

// Une figure s'annonce en posant son ``data-nom`` -- ce que fait de toute
// facon tout trace qui veut exporter un CSV. La barre suit donc les figures
// a venir sans qu'on revienne allonger une liste de selecteurs.
function surveillerFigures() {
  const passe = () => $$("svg.graphe, #carte-modes, svg[data-nom]").forEach(barreExport);
  new MutationObserver(passe).observe(document.body, { childList: true, subtree: true });
  passe();
}

// =====================================================================
// transport / distribution
// =====================================================================
async function calculerImpact() {
  const n = $("#impact-synthese");
  if (!E.reseau) {
    n.textContent = "";
    n.append(el("div", { class: "pb erreur" },
      "Ce reseau n'a pas de zone. L'impact se mesure sur un reseau decrit par zones " +
      "-- transport + instances de distribution -- comme wscc3sm_zones."));
    return;
  }
  const v = E.version;
  $("#impact-etat").textContent = "calcul (quelques secondes par poste)...";
  try {
    const d = await api("/api/impact", {
      reseau: corpsStructure(),
      f_max: parseFloat($("#impact-fmax").value) || 5,
      par_poste: $("#impact-poste").checked,
    });
    if (v !== E.version) return;
    E.res.impact = d;
    $("#impact-etat").textContent = "";
    rendreImpact();
  } catch (err) {
    $("#impact-etat").textContent = "";
    n.textContent = "";
    n.append(el("div", { class: "pb erreur" }, err.message));
  }
}

// --- transport et distribution : GFL au hasard, modeles reduits ------------
function exigerZones(n) {
  if (E.reseau) return true;
  n.textContent = "";
  n.append(el("div", { class: "pb erreur" }, "Il faut un reseau a zones, par exemple wscc9_cigre6 ou ieee14_ieee33."));
  return false;
}

async function ajouterGflAleatoires() {
  const n = $("#mul-ajouts");
  if (!exigerZones(n)) return;
  const g = $("#mul-graine").value;
  try {
    const d = await api("/api/multi/ajouter_gfl", {
      reseau: corpsStructure(),
      n: parseInt($("#mul-n").value, 10) || 1,
      part: (parseFloat($("#mul-part").value) || 5) / 100,
      graine: g === "" ? null : parseInt(g, 10),
    });
    const vue = E.vue;
    entrerZones(d);
    if (vue && vue !== "transport" && E.reseau.zones[vue]) ouvrirVue(vue);
    E.gflAjoutes = [...(E.gflAjoutes || []), ...d.ajouts];
    $("#mul-gfl-etat").textContent = `${d.ajouts.length} GFL ajoute(s), ${E.gflAjoutes.length} depuis le debut`;
    tableau("#mul-ajouts", ["zone", "bus (local)", "P [MW]"],
      E.gflAjoutes.map((a) => ({ cellules: [a.zone, String(a.bus), a.P_W / 1e6] })));
  } catch (err) {
    n.textContent = "";
    n.append(el("div", { class: "pb erreur" }, err.message));
  }
}

async function comparerModelesReduits() {
  const n = $("#mul-synthese");
  if (!exigerZones(n)) return;
  $("#mul-comparer-etat").textContent = "calcul (quelques secondes par distribution)...";
  try {
    const d = await api("/api/multi/comparer", { reseau: corpsStructure() });
    $("#mul-comparer-etat").textContent =
      `production GFL de la distribution : ${fmt(d.p_gfl / 1e6, 1)} MW ; marge du modele complet : max Re = ${fmt(d.max_re, 4)}`;
    const pire = (c) => (c.pire ? `${fmt(c.pire.f, 3)} Hz / ${fmt(c.pire.zeta, 2)} %  ->  ` +
      (c.pire.f_reduit === null ? "absent" : `${fmt(c.pire.f_reduit, 3)} Hz / ${fmt(c.pire.zeta_reduit, 2)} %`) : "");
    const lignes_ = [
      { c: d.transport_rl, nom: "Transport seul, distributions en charges R-L" },
      { c: d.transport_hybride, nom: "Transport seul, equivalent hybride (GFL en puissance constante)" },
      ...d.distributions.map((x) => ({ c: x, nom: `Distribution ${x.instance} seule (bus infini + Zth)` })),
      // la partition de dx/dt = Ax + Bu : l'autre cote passe dans les entrees
      { c: d.partition.transport_quasi_statique, nom: "Partition : transport, x_D quasi statique (A_TT - A_TD A_DD^-1 A_DT)" },
      { c: d.partition.transport_fige, nom: "Partition : transport, x_D fige dans les entrees (A_TT)" },
      { c: d.partition.distribution_fige, nom: "Partition : distributions, x_T fige dans les entrees (A_DD)" },
      { c: d.partition.distribution_quasi_statique, nom: "Partition : distributions, x_T quasi statique" },
    ];
    tableau("#mul-synthese", ["modele reduit", "ecart max [%]", "d zeta max [pts]", "modes absents", "le plus deplace (complet -> reduit)"],
      lignes_.map((l) => ({ l, cellules: [l.nom, 100 * l.c.ecart_max, l.c.dzeta_max, String(l.c.manquants), pire(l.c)] })),
      (x) => (x.l.c.ecart_max > 0.05 || x.l.c.dzeta_max > 2 ? "voisin" : ""));
    const inter = d.interaction.filter((m) => m.zeta < 30);
    tableau("#mul-interaction", ["f [Hz]", "zeta [%]", "part du transport [%]", "elements"],
      inter.map((m) => ({ cellules: [m.f, m.zeta, 100 * m.part_transport, m.elements.map((e) => `${e.nom} ${Math.round(100 * e.part)}%`).join(", ")] })));
  } catch (err) {
    $("#mul-comparer-etat").textContent = "";
    n.textContent = "";
    n.append(el("div", { class: "pb erreur" }, err.message));
  }
}

// Trois series, une couleur et une forme chacune, etiquetees en bout de ligne.
const SERIES_PEN = [
  { cle: "transport_rl", nom: "transport seul (charges R-L)", couleur: "#2a78d6", forme: "rond" },
  { cle: "transport_hybride", nom: "transport seul (hybride)", couleur: "#eb6834", forme: "carre" },
  { cle: "distribution", nom: "distribution seule (IB + Zth)", couleur: "#1baf7a", forme: "losange" },
];

async function balayerPenetration() {
  const n = $("#pen-table");
  if (!exigerZones(n)) return;
  $("#pen-etat").textContent = "calcul (une dizaine de secondes par palier)...";
  try {
    const d = await api("/api/multi/penetration", {
      reseau: corpsStructure(),
      n_max: parseInt($("#pen-nmax").value, 10) || 40,
      pas: parseInt($("#pen-pas").value, 10) || 10,
      graine: parseInt($("#pen-graine").value, 10) || 1,
    });
    E.res.penetration = d;
    const ok = d.paliers.filter((p) => !p.erreur);
    const ko = d.paliers.find((p) => p.erreur);
    $("#pen-etat").textContent = `${ok.length} palier(s)` + (ko ? ` ; arret a ${ko.n_gfl} GFL : ${ko.erreur}` : "");
    dessinerPenetration();
    tableau("#pen-table", ["GFL ajoutes", "production GFL [MW]", "max Re", ...SERIES_PEN.flatMap((x) => [`${x.nom} : ecart [%]`, "d zeta [pts]"]), "modes d'interaction"],
      ok.map((p) => ({ cellules: [String(p.n_gfl), p.p_gfl / 1e6, p.max_re, ...SERIES_PEN.flatMap((x) => [100 * p[x.cle].ecart, p[x.cle].dzeta]), String(p.interaction)] })));
  } catch (err) {
    $("#pen-etat").textContent = "";
    n.textContent = "";
    n.append(el("div", { class: "pb erreur" }, err.message));
  }
}

function dessinerPenetration() {
  const d = E.res.penetration;
  const svg = $("#carte-penetration");
  svg.textContent = "";
  const ok = d ? d.paliers.filter((p) => !p.erreur) : [];
  if (!ok.length) return;
  const xs = ok.map((p) => p.p_gfl / 1e6);
  const ys = ok.flatMap((p) => SERIES_PEN.map((x) => 100 * p[x.cle].ecart));
  const [x0, x1] = minMax(xs), [, y1] = minMax(ys);
  const b = { x0: x0 - 0.02 * (x1 - x0 || 1), x1: x1 + 0.25 * (x1 - x0 || 1), y0: 0, y1: (y1 || 1) * 1.12 };
  const P = planLineaire(svg, b, { W: 1000, H: 400, ML: 70, MR: 18, MT: 30, MB: 44, titreX: "production GFL de la distribution [MW]", titreY: "ecart relatif max |d lambda| / |lambda| [%]" });
  const marque = (x, y, s, attrs) => s.forme === "rond"
    ? el("circle", { cx: x, cy: y, r: 5, ...attrs })
    : s.forme === "carre"
      ? el("rect", { x: x - 4.5, y: y - 4.5, width: 9, height: 9, rx: 1.5, ...attrs })
      : el("rect", { x: x - 4.5, y: y - 4.5, width: 9, height: 9, transform: `rotate(45 ${x} ${y})`, ...attrs });
  for (const s of SERIES_PEN) {
    const pts = ok.map((p) => [P.PX(p.p_gfl / 1e6), P.PY(100 * p[s.cle].ecart)]);
    P.g.append(el("polyline", { points: pts.map((q) => q.join(",")).join(" "), fill: "none", stroke: s.couleur, "stroke-width": 2 }));
    ok.forEach((p, i) => {
      const m = marque(pts[i][0], pts[i][1], s, { fill: s.couleur, stroke: "#fff", "stroke-width": 1.5 });
      infobulle(m, () => `${s.nom}\n${p.n_gfl} GFL ajoutes, ${fmt(p.p_gfl / 1e6, 1)} MW\necart ${fmt(100 * p[s.cle].ecart, 3)} %, d zeta ${fmt(p[s.cle].dzeta, 2)} pts`);
      P.g.append(m);
    });
    const [lx, ly] = pts[pts.length - 1];
    P.g.append(el("text", { x: lx + 10, y: ly + 4, class: "etiquette legende-svg" }, s.nom));
  }
  svg.dataset.nom = "penetration";
  svg.dataset.csv = "n_gfl,p_gfl_MW,max_re," + SERIES_PEN.flatMap((x) => [`${x.cle}_ecart_pct`, `${x.cle}_dzeta`]).join(",") + "\n" +
    ok.map((p) => [p.n_gfl, p.p_gfl / 1e6, p.max_re, ...SERIES_PEN.flatMap((x) => [100 * p[x.cle].ecart, p[x.cle].dzeta])].join(",")).join("\n");
}

function tableau(cible, entetes, lignes_, classe) {
  const n = $(cible);
  n.textContent = "";
  if (!lignes_.length) {
    n.append(el("p", { class: "doux" }, "aucun"));
    return;
  }
  const t = el("table");
  const th = el("tr");
  entetes.forEach((x) => th.append(el("th", {}, x)));
  t.append(el("thead", {}, th));
  const tb = el("tbody");
  for (const l of lignes_) {
    const tr = el("tr", { class: classe ? classe(l) : "" });
    for (const v of l.cellules) tr.append(el("td", { class: typeof v === "string" ? "mono" : "chiffre" }, typeof v === "string" ? v : fmt(v, 3)));
    tb.append(tr);
  }
  t.append(tb);
  n.append(t);
}

function rendreImpact() {
  const d = E.res.impact;
  if (!d) return;
  const n = $("#impact-synthese");
  n.textContent = "";
  const Sb = E.plat ? E.plat.tables.Sb : E.tables.Sb;
  const maxMove = d.paires.reduce((a, m) => Math.max(a, Math.abs(m.ecart_zeta)), 0);
  const phrases = [
    `Point de fonctionnement du transport identique a ${d.ecart_tension.toExponential(1)} pu pres : seule la dynamique differe.`,
    `Marge : max Re = ${fmt(d.max_re_equivalent, 4)} sans la dynamique de la distribution, ${fmt(d.max_re_complet, 4)} avec.`,
    `Plus grand deplacement d'amortissement d'un mode du transport : ${fmt(maxMove, 2)} point(s).`,
    `${d.crees.length} mode(s) cree(s), ${d.degrades.length} mode(s) degrade(s) par la distribution.`,
  ];
  for (const ph of phrases) n.append(el("p", {}, ph));
  const t = el("table");
  t.append(el("thead", {}, el("tr", {}, el("th", {}, "poste"), el("th", {}, "P soutire [MW]"), el("th", {}, "Q [Mvar]"))));
  const tb = el("tbody");
  for (const [k, s2] of Object.entries(d.soutirage))
    tb.append(el("tr", {}, el("td", { class: "mono" }, k), el("td", { class: "chiffre" }, fmt(s2.P * Sb / 1e6, 1)), el("td", { class: "chiffre" }, fmt(s2.Q * Sb / 1e6, 1))));
  t.append(tb);
  n.append(t);

  tableau("#impact-paires",
    ["f sans [Hz]", "zeta sans [%]", "f avec [Hz]", "zeta avec [%]", "ecart f [%]", "ecart zeta [pts]", "part distribution [%]", "etat dominant"],
    d.paires.map((m) => ({ m, cellules: [m.f_eq, m.zeta_eq, m.f, m.zeta, m.ecart_f, m.ecart_zeta, m.part_distribution, m.dominant] })),
    (l) => (Math.abs(l.m.ecart_zeta) > 1 || l.m.part_distribution > 10 ? "voisin" : ""));
  tableau("#impact-alertes",
    ["nature", "f [Hz]", "zeta [%]", "part distribution [%]", "etat dominant"],
    [...d.crees.map((c) => ({ cellules: ["cree", c.f, c.zeta, c.part_distribution, c.dominant] })),
     ...d.degrades.map((m) => ({ cellules: [`degrade (zeta ${fmt(m.zeta_eq, 2)} % sans)`, m.f, m.zeta, m.part_distribution, m.dominant] }))],
    (l) => (l.cellules[2] < 0 ? "instable" : l.cellules[2] < 5 ? "voisin" : ""));
  const postes = Object.entries(d.par_poste);
  if (postes.length) {
    const fs = postes[0][1].map((x) => x.f_eq);
    tableau("#impact-postes",
      ["poste (seul dynamique)", ...fs.map((f) => `${fmt(f, 2)} Hz : df % / dzeta pts`)],
      postes.map(([k, v]) => ({ cellules: [k, ...v.map((x) => `${fmt(x.ecart_f, 2)} / ${fmt(x.ecart_zeta, 2)}`)] })));
  } else $("#impact-postes").textContent = "";

  dessinerImpact();
}

// La carte du transport : creux = equivalent statique, plein = modele complet,
// un trait relie les deux. Meme plan lineaire que les autres cartes des modes.
const IMPACT_LEGENDE = [
  { forme: "creux", couleur: "#1f4e79", texte: "distribution en charge" },
  { forme: "plein", couleur: "#1f4e79", texte: "modele complet" },
  { forme: "plein", couleur: "#7a4fa3", texte: "complet, distribution > 10 %" },
  { forme: "losange", couleur: "#bf2617", texte: "mode cree" },
];

function dessinerImpact() {
  const d = E.res.impact;
  if (!d) return;
  const svg = $("#carte-impact");
  svg.textContent = "";
  barreBornes("imp", dessinerImpact);
  const fmax = parseFloat($("#impact-fmax").value) || 5;
  const crees = d.crees.filter((c) => c.f <= fmax * 1.5);
  const pts = [...d.paires.flatMap((m) => [[m.re_eq, m.f_eq], [m.re, m.f]]), ...crees.map((c) => [c.re, c.f])];
  if (!pts.length) return;
  const b = bornesPlan("imp", pts.map((q) => q[0]), pts.map((q) => q[1]), { zeroX: true, zeroY: true });
  const P = planLineaire(svg, b, { W: 1000, H: 480, ML: 70, MR: 18, MT: 44, MB: 44, titreX: "Re(lambda) [1/s]", titreY: "f [Hz]" });

  // legende, au-dessus du trace
  let lx = 76;
  const marque = (x, y, e) => e.forme === "creux"
    ? el("circle", { cx: x, cy: y, r: 5.5, fill: "#fff", stroke: e.couleur, "stroke-width": 2 })
    : e.forme === "plein"
      ? el("circle", { cx: x, cy: y, r: 5, fill: e.couleur })
      : el("rect", { x: x - 5, y: y - 5, width: 10, height: 10, fill: e.couleur, transform: `rotate(45 ${x} ${y})` });
  for (const e of IMPACT_LEGENDE) {
    svg.append(marque(lx, 12, e));
    const t = el("text", { x: lx + 10, y: 16, class: "etiquette legende-svg" }, e.texte);
    svg.append(t);
    lx += 12 + e.texte.length * 6.2 + 26;
  }
  svg.append(el("line", { x1: lx, y1: 12, x2: lx + 22, y2: 12, stroke: "#98a2b3", "stroke-width": 1.5 }));
  svg.append(el("text", { x: lx + 28, y: 16, class: "etiquette legende-svg" }, "deplacement du mode"));

  for (const m of d.paires) {
    P.g.append(el("line", { x1: P.PX(m.re_eq), y1: P.PY(m.f_eq), x2: P.PX(m.re), y2: P.PY(m.f), stroke: "#98a2b3", "stroke-width": 1.5 }));
    P.g.append(marque(P.PX(m.re_eq), P.PY(m.f_eq), IMPACT_LEGENDE[0]));
    const c = marque(P.PX(m.re), P.PY(m.f), IMPACT_LEGENDE[m.part_distribution > 10 ? 2 : 1]);
    infobulle(c, () => `${fmt(m.f_eq, 3)} Hz / ${fmt(m.zeta_eq, 2)} %  ->  ${fmt(m.f, 3)} Hz / ${fmt(m.zeta, 2)} %` + "\n" +
      `distribution : ${fmt(m.part_distribution, 1)} % ; ${m.dominant}`);
    P.g.append(c);
  }
  for (const c of crees) {
    const q = marque(P.PX(c.re), P.PY(c.f), IMPACT_LEGENDE[3]);
    infobulle(q, () => `CREE : ${fmt(c.f, 3)} Hz, zeta ${fmt(c.zeta, 2)} %` + "\n" + `distribution ${fmt(c.part_distribution, 1)} % ; ${c.dominant}`);
    P.g.append(q);
  }
  zoomPlan(svg, "imp", P);
  svg.dataset.nom = "impact_distribution";
  svg.dataset.csv = "f_sans,zeta_sans,f_avec,zeta_avec,part_distribution,dominant" + "\n" +
    d.paires.map((m) => [m.f_eq, m.zeta_eq, m.f, m.zeta, m.part_distribution, m.dominant].join(",")).join("\n");
}

// =====================================================================
// balayage
// =====================================================================
// Un balayage repond a la question que ni une carte des modes ni une
// simulation ne posent : a partir de QUELLE valeur ce reglage rend le
// reseau instable. D'ou les deux traces -- le lieu des racines, qui dit par
// ou les modes passent, et la marge, qui dit quand.

// Le formulaire n'a de sens que par morceaux : un reglage de famille n'a pas
// de bus, une cellule de table n'a pas de reglage. On montre ce qui sert.
// Le parametre a balayer se choisit comme on le pense : un ELEMENT (SM2,
// L3_5_7, la charge du bus 5...), puis l'un de ses parametres. Chaque
// parametre porte sa cible pour le serveur et sa valeur actuelle.
const COLONNES_TOPOLOGIE = new Set(["ID", "BusID", "From", "To", "BusAC", "BusDC", "Type", "Mode"]);
const PREFIXE_UNITE = { T_SMs: "SM", T_GFM1s: "GFM", T_GFM2s: "GFM", T_GFM3s: "GFM", T_GFM4s: "GFM", T_GFL1s: "GFL", T_GFL2s: "GFL", T_IBs: "IB" };

/** Les tables sur lesquelles le serveur calcule (le reseau aplati s'il a des zones). */
function tablesCalcul() {
  return E.reseau && E.plat ? E.plat.tables : E.tables;
}

function elementsBalayage() {
  const T = tablesCalcul();
  const rows = (t) => T[t] || [];
  const out = [];
  const colonnes = (table, i, r) =>
    COLS(table).map((c, k) => ({ c, k })).filter(({ c, k }) => !COLONNES_TOPOLOGIE.has(c) && typeof r[k] === "number")
      .map(({ c, k }) => ({ nom: c, valeur: r[k], cible: { genre: "table", table, ligne: i, colonne: c } }));
  const yder = (bus) => rows("Y_DER").findIndex((r) => +r[0] === bus);
  const reglages = (fam, bus, table) => {
    const tab = { SM: "SM_PARAMS", GFM: "GFM_PARAMS", GFL: "GFL_PARAMS" }[fam];
    const eff = Object.assign({}, defautsReglages(fam, bus, table), T[tab] || {}, (T.PARAMS_PAR_BUS || {})[String(bus)] || {});
    return Object.entries(eff).filter(([, v]) => typeof v === "number")
      .map(([k, v]) => ({ nom: k, valeur: v, cible: { genre: "bus", bus, cle: k } }));
  };
  // unites
  for (const [table, pre] of Object.entries(PREFIXE_UNITE)) {
    rows(table).forEach((r, i) => {
      const bus = +r[1];
      const params = colonnes(table, i, r);
      const iy = yder(bus);
      if (iy >= 0) {
        const y = rows("Y_DER")[iy];
        params.push(
          { nom: "P consigne (Pg_W)", valeur: y[5], cible: { genre: "table", table: "Y_DER", ligne: iy, colonne: "Pg_W" } },
          { nom: "Q (Qg_var)", valeur: y[6], cible: { genre: "table", table: "Y_DER", ligne: iy, colonne: "Qg_var" } },
          { nom: "V consigne (V_pu)", valeur: y[3], cible: { genre: "table", table: "Y_DER", ligne: iy, colonne: "V_pu" } });
      }
      const fam = familleReglages(table);
      if (fam) params.push(...reglages(fam, bus, table));
      out.push({ cle: `${table}:${i}`, nom: `${pre}${bus}`, params });
    });
  }
  // branches
  for (const table of ["T_Lines_pu", "T_Lines_SI"])
    rows(table).forEach((r, i) => out.push({ cle: `${table}:${i}`, nom: `L${r[0]}_${r[1]}_${r[2]}`, params: colonnes(table, i, r) }));
  rows("T_Transfos").forEach((r, i) => out.push({ cle: `T_Transfos:${i}`, nom: `T${r[0]}_${r[1]}_${r[2]}`, params: colonnes("T_Transfos", i, r) }));
  // charges et shunts
  rows("Y_DER").forEach((r, i) => {
    if (!(+r[7] || +r[8])) return;
    out.push({
      cle: `Y_DER:${i}`, nom: `Charge du bus ${r[0]}`, params: [
        { nom: "P (Pl_W)", valeur: r[7], cible: { genre: "table", table: "Y_DER", ligne: i, colonne: "Pl_W" } },
        { nom: "Q (Ql_var)", valeur: r[8], cible: { genre: "table", table: "Y_DER", ligne: i, colonne: "Ql_var" } },
      ],
    });
  });
  rows("T_Loads").forEach((r, i) => out.push({ cle: `T_Loads:${i}`, nom: `Load${r[4]}${r[0]} (bus ${r[1]})`, params: colonnes("T_Loads", i, r) }));
  rows("T_Shunts").forEach((r, i) => out.push({ cle: `T_Shunts:${i}`, nom: `${+r[2] < 0 ? "Inductance" : "Condensateur"} ${r[0]} (bus ${r[1]})`, params: colonnes("T_Shunts", i, r) }));
  // continu
  rows("T_BusDC").forEach((r, i) => out.push({ cle: `T_BusDC:${i}`, nom: `Bus DC ${r[0]}`, params: colonnes("T_BusDC", i, r) }));
  rows("T_LinesDC").forEach((r, i) => out.push({ cle: `T_LinesDC:${i}`, nom: `LDC${r[0]}_${r[1]}_${r[2]}`, params: colonnes("T_LinesDC", i, r) }));
  rows("T_ShuntsDC").forEach((r, i) => out.push({ cle: `T_ShuntsDC:${i}`, nom: `Condensateur DC ${r[0]} (bus ${r[1]})`, params: colonnes("T_ShuntsDC", i, r) }));
  rows("T_DCDC").forEach((r, i) => out.push({ cle: `T_DCDC:${i}`, nom: `DCDC${r[0]} ${r[3]} (bus ${r[1]})`, params: colonnes("T_DCDC", i, r) }));
  rows("T_LoadsDC").forEach((r, i) => out.push({ cle: `T_LoadsDC:${i}`, nom: `Charge DC ${r[0]} (bus ${r[1]})`, params: colonnes("T_LoadsDC", i, r) }));
  // reglages communs a toute une famille
  for (const [fam, tab] of [["SM", "SM_PARAMS"], ["GFM", "GFM_PARAMS"], ["GFL", "GFL_PARAMS"]]) {
    const eff = Object.assign({}, E.defauts.reglages[fam] || {}, T[tab] || {});
    out.push({
      cle: `famille:${fam}`, nom: `Tous les ${fam} (reglage de famille)`,
      params: Object.entries(eff).filter(([, v]) => typeof v === "number")
        .map(([k, v]) => ({ nom: k, valeur: v, cible: { genre: "famille", famille: tab, cle: k } })),
    });
  }
  return out.filter((e) => e.params.length);
}

/** Remplit la liste des elements, en gardant le choix courant s'il existe. */
function rendreCibleBalayage() {
  if (!E.defauts || !E.tables) return;
  E.elementsBal = elementsBalayage();
  const s = $("#bal-element");
  const garde = s.value;
  s.textContent = "";
  for (const e of E.elementsBal) s.append(el("option", { value: e.cle }, e.nom));
  if (E.elementsBal.some((e) => e.cle === garde)) s.value = garde;
  remplirParametres(false);
}

function elementBalayage() {
  return (E.elementsBal || []).find((e) => e.cle === $("#bal-element").value);
}

function parametreBalayage() {
  const e = elementBalayage();
  return e ? e.params[+$("#bal-param").value] : null;
}

/** Les parametres de l'element. ``nouveau`` : l'element vient de changer. */
function remplirParametres(nouveau) {
  const e = elementBalayage();
  const s = $("#bal-param");
  const garde = s.selectedOptions[0] ? s.selectedOptions[0].textContent : "";
  s.textContent = "";
  if (!e) return;
  e.params.forEach((q, i) => s.append(el("option", { value: i }, q.nom)));
  const i = e.params.findIndex((q) => q.nom === garde);
  const change = nouveau || i < 0 || E.parametreBal !== `${e.cle}/${garde}`;
  s.value = i >= 0 && !nouveau ? i : 0;
  if (change) choisirParametre();
  else afficherValeurActuelle();
}

function afficherValeurActuelle() {
  const q = parametreBalayage();
  $("#bal-actuelle").textContent = q ? fmt(q.valeur, 6) : "";
}

/** Un NOUVEAU parametre : De et A reprennent leurs valeurs par defaut,
 *  autour de la valeur actuelle (de la moitie a une fois et demie). */
function choisirParametre() {
  const q = parametreBalayage(), e = elementBalayage();
  afficherValeurActuelle();
  if (!q) return;
  E.parametreBal = `${e.cle}/${q.nom}`;
  const v = +q.valeur;
  let [a, b] = v ? [0.5 * v, 1.5 * v] : [0, 1];
  if (a > b) [a, b] = [b, a];
  $("#bal-min").value = String(+a.toPrecision(4));
  $("#bal-max").value = String(+b.toPrecision(4));
}

function cibleBalayage() {
  const q = parametreBalayage();
  if (!q) throw new Error("choisir un element et un parametre");
  return q.cible;
}

// ------------------------------------------------ robustesse des modes
async function robustesse() {
  const v = E.version;
  $("#rob-etat").textContent = "calcul des variantes...";
  for (const c of ["#rob-alertes", "#rob-modes", "#rob-detail"]) $(c).textContent = "";
  try {
    const d = await lancerTache("robustesse", corpsReseau({
      frac: (parseFloat($("#rob-frac").value) || 10) / 100,
      // TOUS les modes : pas de bande, l'utilisateur zoome lui-meme
      f_min: 0,
      f_max: 1e12,
      zeta_max: null,
      puissance: $("#rob-p").checked,
      deconnexion: $("#rob-n1").checked,
    }));
    if (v !== E.version) return;
    E.res.robustesse = d;
    E.robSel = 0;
    $("#rob-etat").textContent = `${d.n_variantes} variantes, ${d.modes.length} mode(s)`;
    rendreRobustesse();
  } catch (err) {
    $("#rob-etat").textContent = "";
    $("#rob-alertes").append(el("div", { class: "pb erreur" }, err.message));
  }
}

function rendreRobustesse() {
  const d = E.res.robustesse;
  if (!d) return;
  const al = $("#rob-alertes");
  al.textContent = "";
  for (const x of d.instables) {
    const qui = x.elements.map((e) => `${e.nom} ${Math.round(100 * e.part)} %`).join(", ");
    al.append(el("div", { class: "pb erreur" }, `Instable apres ${x.libelle} : max Re = ${fmt(x.max_re, 3)} (${qui})`));
  }
  for (const x of d.echecs)
    al.append(el("div", { class: "pb avis" },
      `Sans point de fonctionnement apres ${x.libelle}${x.genre === "ilotage" ? " (ilotage)" : ""} : ${x.raison}`));
  const n = $("#rob-modes");
  n.textContent = "";
  if (!d.modes.length) {
    n.append(el("p", { class: "doux" }, "aucun mode"));
    $("#carte-robustesse").textContent = "";
    $("#rob-detail").textContent = "";
    return;
  }
  const t = el("table");
  t.append(el("thead", {}, el("tr", {},
    ...["f [Hz]", "zeta [%]", "elements", "score P", "score N-1", "zeta min [%]", "pire perturbation", "disparait avec"]
      .map((h) => el("th", {}, h)))));
  const tb = el("tbody");
  d.modes.forEach((m, i) => {
    tb.append(el("tr", {
      class: [m.zeta_min < 0 ? "instable" : m.zeta_min < 5 ? "voisin" : "", i === E.robSel ? "choisi" : ""].join(" "),
      style: "cursor:pointer",
      onclick: () => choisirModeRobustesse(i),
    },
      el("td", { class: "chiffre" }, fmt(m.f, 3)),
      el("td", { class: "chiffre" }, fmt(m.zeta, 3)),
      el("td", { class: "mono" }, m.elements.map((e) => `${e.nom} ${Math.round(100 * e.part)}%`).join(", ")),
      el("td", { class: "chiffre" }, fmt(m.score_puissance, 3)),
      el("td", { class: "chiffre" }, fmt(m.score_deconnexion, 3)),
      el("td", { class: m.zeta_min < 0 ? "chiffre instable" : "chiffre" }, fmt(m.zeta_min, 3)),
      el("td", { class: "mono" }, m.pire ? m.pire.libelle : ""),
      el("td", { class: "mono" }, m.disparitions.join(", "))));
  });
  t.append(tb);
  n.append(t);
  dessinerRobustesse();
  detailRobustesse();
}

function choisirModeRobustesse(i) {
  E.robSel = i;
  rendreRobustesse();
}

// Couleurs par GENRE de perturbation (identite, ordre fixe), et une forme en
// plus : rond = puissance, carre = deconnexion. L'identite ne tient jamais a la
// seule couleur.
const ROB_GENRES = {
  puissance: { couleur: "#2a78d6", nom: "variation de puissance" },
  deconnexion: { couleur: "#eb6834", nom: "deconnexion (N-1)" },
};

function dessinerRobustesse() {
  const d = E.res.robustesse;
  if (!d) return;
  const svg = $("#carte-robustesse");
  svg.textContent = "";
  barreBornes("rob", dessinerRobustesse);
  // le plan complexe entier : chaque mode et son conjugue
  const garde = (x) => !x.disparu && x.genre !== "ilotage";
  const xs = [], ys = [];
  for (const m of d.modes) {
    xs.push(m.re); ys.push(m.im, -m.im);
    for (const x of m.effets) if (garde(x)) { xs.push(x.re); ys.push(x.im, -x.im); }
  }
  const b = bornesPlan("rob", xs, ys, { zeroX: true, zeroY: true });
  const P = planLineaire(svg, b, { W: 1000, H: 480, ML: 70, MR: 18, MT: 30, MB: 44, titreX: "Re(lambda) [1/s]", titreY: "Im(lambda) [rad/s]" });
  isoAmortissement(P, false, [0.05, 0.1]);
  const dedans = (re, im) => re >= b.x0 && re <= b.x1 && im >= b.y0 && im <= b.y1;

  // legende, au-dessus du trace
  let lx = 1000 - 18 - 470;
  for (const [genre, st] of Object.entries(ROB_GENRES)) {
    svg.append(genre === "puissance"
      ? el("circle", { cx: lx, cy: 12, r: 5, fill: st.couleur })
      : el("rect", { x: lx - 5, y: 7, width: 10, height: 10, rx: 1.5, fill: st.couleur }));
    svg.append(el("text", { x: lx + 10, y: 16, class: "etiquette" }, st.nom));
    lx += 170;
  }
  svg.append(el("circle", { cx: lx, cy: 12, r: 6, fill: "#1d2433", stroke: "#fff", "stroke-width": 2 }));
  svg.append(el("text", { x: lx + 10, y: 16, class: "etiquette" }, "cas de base"));

  // les variantes de tous les modes, en un chemin par genre (des milliers de
  // points : un element SVG chacun serait trop lourd), estompees
  const sel = E.robSel;
  const chemins = { puissance: "", deconnexion: "" };
  d.modes.forEach((m, i) => {
    if (i === sel) return;
    for (const x of m.effets) {
      if (!garde(x)) continue;
      for (const im of [x.im, -x.im]) {
        if (!dedans(x.re, im)) continue;
        const px = P.PX(x.re).toFixed(1), py = P.PY(im).toFixed(1);
        chemins[x.genre] += x.genre === "puissance"
          ? `M${+px + 4},${py}a4,4 0 1,0 -8,0a4,4 0 1,0 8,0`
          : `M${+px - 4},${+py - 4}h8v8h-8z`;
      }
    }
  });
  for (const [genre, dd] of Object.entries(chemins))
    if (dd) P.g.append(el("path", { d: dd, fill: ROB_GENRES[genre].couleur, "fill-opacity": 0.25 }));

  const marque = (x, y, genre, attrs) => genre === "puissance"
    ? el("circle", { cx: x, cy: y, r: 5, ...attrs })
    : el("rect", { x: x - 5, y: y - 5, width: 10, height: 10, rx: 1.5, ...attrs });

  // les modes du cas de base (et leurs conjugues)
  d.modes.forEach((m, i) => {
    for (const im of m.im > 1e-9 ? [m.im, -m.im] : [m.im]) {
      if (!dedans(m.re, im)) continue;
      const n = el("circle", { cx: P.PX(m.re), cy: P.PY(im), r: i === sel ? 7 : 5, fill: "#1d2433", stroke: "#fff", "stroke-width": 1.5, style: "cursor:pointer" });
      infobulle(n, () => `mode ${fmt(m.re, 4)} ${m.im >= 0 ? "+" : "-"} ${fmt(Math.abs(m.im), 4)}j\n${fmt(m.f, 3)} Hz, zeta ${fmt(m.zeta, 2)} %\n` +
        `${m.elements.map((e) => `${e.nom} ${Math.round(100 * e.part)} %`).join(", ")}\nzeta min ${fmt(m.zeta_min, 2)} %${m.pire ? " (" + m.pire.libelle + ")" : ""}`);
      n.addEventListener("click", () => choisirModeRobustesse(i));
      P.g.append(n);
    }
  });

  // le mode choisi : ses variantes, reliees a lui, au premier plan
  const m = d.modes[sel];
  if (m) {
    for (const x of m.effets) {
      if (!garde(x)) continue;
      const col = ROB_GENRES[x.genre].couleur;
      for (const [imb, im] of [[m.im, x.im], [-m.im, -x.im]]) {
        if (imb === -m.im && m.im <= 1e-9) continue;
        P.g.append(el("line", { x1: P.PX(m.re), y1: P.PY(imb), x2: P.PX(x.re), y2: P.PY(im), stroke: col, "stroke-width": 1, opacity: 0.5 }));
        const n = marque(P.PX(x.re), P.PY(im), x.genre, { fill: col, stroke: "#fff", "stroke-width": 1.5 });
        infobulle(n, () => `${x.libelle}\n${fmt(x.re, 4)} ${x.im >= 0 ? "+" : "-"} ${fmt(Math.abs(x.im), 4)}j\n` +
          `f = ${fmt(x.f, 3)} Hz (${x.df >= 0 ? "+" : ""}${fmt(x.df, 3)})\nzeta = ${fmt(x.zeta, 2)} % (${x.dzeta >= 0 ? "+" : ""}${fmt(x.dzeta, 2)} pts)`);
        P.g.append(n);
      }
    }
    for (const im of m.im > 1e-9 ? [m.im, -m.im] : [m.im])
      P.g.append(el("circle", { cx: P.PX(m.re), cy: P.PY(im), r: 7, fill: "#1d2433", stroke: "#fff", "stroke-width": 2, "pointer-events": "none" }));
  }
  zoomPlan(svg, "rob", P);
  svg.dataset.nom = "robustesse";
  svg.dataset.csv = "mode_Re,mode_Im,mode_f_Hz,mode_zeta,perturbation,genre,Re,Im,f_Hz,zeta,df_Hz,dzeta_pts,score\n" +
    d.modes.flatMap((mm) => mm.effets.map((x) => [mm.re, mm.im, mm.f, mm.zeta, `"${x.libelle}"`, x.genre, x.re, x.im, x.f, x.zeta, x.df, x.dzeta, x.score].join(","))).join("\n");
}

function detailRobustesse() {
  const m = E.res.robustesse.modes[E.robSel];
  $("#rob-titre").textContent = m ? `${fmt(m.f, 3)} Hz, zeta ${fmt(m.zeta, 2)} %` : "";
  if (!m) return;
  tableau("#rob-detail",
    ["perturbation", "genre", "f [Hz]", "zeta [%]", "df [Hz]", "dzeta [pts]", "score"],
    m.effets.map((x) => ({
      x,
      cellules: x.disparu
        ? [x.libelle, x.genre, "disparu", "", "", "", ""]
        : [x.libelle, x.genre, x.f, x.zeta, x.df, x.dzeta, x.score],
    })),
    (l) => (l.x.disparu || l.x.genre === "ilotage" ? "doux" : l.x.zeta < 0 ? "instable" : l.x.dzeta < -2 ? "voisin" : ""));
}

// ------------------------------------------- hypotheses de simplification
// E.simp : {type ou equipement: {etat: "statique" | "fige" | "dynamique"}},
// la seule source de verite de l'onglet ; les preregles ne font qu'y ecrire.
const SIMP_SORTS = [["", "dynamique"], ["statique", "dx/dt = 0"], ["fige", "fige"]];

async function chargerSimplification() {
  const v = E.version;
  $("#simp-etat").textContent = "lecture du reseau...";
  try {
    const d = await api("/api/simplification/inventaire", corpsReseau());
    if (v !== E.version) return;
    E.simpInv = d;
    E.simpInvVersion = v;
    // les reglages d'un autre reseau ne valent pas pour celui-ci
    const connus = new Set(d.types.flatMap((t) => [t.type, ...t.elements.map((e) => e.nom)]));
    for (const k of Object.keys(E.simp)) if (!connus.has(k)) delete E.simp[k];
    // un reseau enregistre avec ses hypotheses les rapporte
    if (!Object.keys(E.simp).length && d.simplifications && Object.keys(d.simplifications).length) {
      E.simp = JSON.parse(JSON.stringify(d.simplifications));
      E.simpAppliquer = true;
    }
    $("#simp-appliquer").checked = !!E.simpAppliquer;
    $("#simp-etat").textContent = `${d.types.length} type(s) d'equipement`;
    rendreSimplificationReglages();
  } catch (err) {
    $("#simp-etat").textContent = "";
    $("#simp-alertes").textContent = "";
    $("#simp-alertes").append(el("div", { class: "pb erreur" }, err.message));
  }
}

function simpNettoyer() {
  for (const [k, v] of Object.entries(E.simp)) {
    for (const [e, s] of Object.entries(v)) if (!s || s === "dynamique") delete v[e];
    if (!Object.keys(v).length) delete E.simp[k];
  }
}

function simpChange() {
  simpNettoyer();
  // les autres calculs dependent du modele : les resultats sont perimes
  if (E.simpAppliquer) {
    E.version += 1;
    E.res = { lf: null, modal: null, signaux: null, forme: null, balayage: null, validation: null, impact: null };
    E.simpInvVersion = E.version;
  }
  E.res.simplification = null;
  rendreSimplificationReglages();
  sauverBrouillon();
}

function simpSelect(cle, etat, vide) {
  const cur = (E.simp[cle] || {})[etat] || "";
  const s = el("select", { class: cur ? "simp-actif" : "" },
    ...SIMP_SORTS.map(([val, lab]) => el("option", { value: val, ...(val === cur ? { selected: "" } : {}) }, val === "" && vide ? vide : lab)));
  s.addEventListener("change", () => {
    (E.simp[cle] = E.simp[cle] || {})[etat] = s.value;
    simpChange();
  });
  return s;
}

function rendreSimplificationReglages() {
  const d = E.simpInv;
  if (!d) return;
  // les preregles : coche si toutes leurs entrees sont en place
  const pr = $("#simp-preregles");
  pr.textContent = "";
  if (!d.preregles.length) pr.append(el("p", { class: "doux" }, "aucun pour ce reseau"));
  for (const p of d.preregles) {
    const tout = Object.entries(p.simplifications).every(([k, v]) =>
      Object.entries(v).every(([e, s]) => (E.simp[k] || {})[e] === s));
    const c = el("input", { type: "checkbox", ...(tout ? { checked: "" } : {}) });
    c.addEventListener("change", () => {
      for (const [k, v] of Object.entries(p.simplifications))
        for (const [e, s] of Object.entries(v)) (E.simp[k] = E.simp[k] || {})[e] = c.checked ? s : "";
      simpChange();
    });
    const types = Object.keys(p.simplifications).join(", ");
    pr.append(el("label", { class: "case-service simp-preregle", title: JSON.stringify(p.simplifications) },
      c, el("span", {}, el("b", {}, p.nom.replace(/_/g, " ")), " : ", p.description, el("span", { class: "doux" }, ` (${types})`))));
  }
  // les types, et sous chacun ses equipements
  const n = $("#simp-types");
  const ouverts = new Set($$("#simp-types details[open]").map((x) => x.dataset.type));
  n.textContent = "";
  for (const t of d.types) {
    const nSimp = Object.keys(E.simp[t.type] || {}).length
      + t.elements.filter((e) => E.simp[e.nom]).length;
    const det = el("details", { "data-type": t.type, ...(ouverts.has(t.type) ? { open: "" } : {}) });
    const etatsType = el("div", { class: "simp-etats" },
      ...t.etats.map((e) => el("label", {}, el("span", { class: "mono" }, e), simpSelect(t.type, e))));
    det.append(el("summary", {}, el("b", {}, t.type), ` ${t.elements.length} equipement(s), ${t.etats.length} etat(s)`,
      nSimp ? el("span", { class: "simp-compte" }, ` ${nSimp} reglage(s)`) : ""));
    det.append(etatsType);
    // les equipements : rendus a l'ouverture seulement (des centaines de lignes)
    const liste = el("div", { class: "simp-elements" });
    const remplir = () => {
      if (liste.childElementCount) return;
      for (const x of t.elements) {
        const famille = Object.entries(x.dynamiques).filter(([, s]) => s !== "dynamique").map(([e, s]) => `${e} ${s}`);
        liste.append(el("div", { class: "simp-element" },
          el("div", { class: "mono" }, el("b", {}, x.nom), famille.length ? el("span", { class: "doux" }, ` (famille : ${famille.join(", ")})`) : ""),
          el("div", { class: "simp-etats" }, ...t.etats.map((e) => el("label", {}, el("span", { class: "mono" }, e), simpSelect(x.nom, e, "comme le type"))))));
      }
    };
    det.addEventListener("toggle", () => det.open && remplir());
    if (det.open) remplir();
    det.append(liste);
    n.append(det);
  }
}

async function comparerSimplification() {
  const v = E.version;
  $("#simp-etat").textContent = "calcul des deux spectres...";
  $("#simp-alertes").textContent = "";
  try {
    simpNettoyer();
    const d = await api("/api/simplification/comparer", corpsReseau({ simplifications: E.simp }));
    if (v !== E.version) return;
    E.res.simplification = d;
    $("#simp-etat").textContent = `${d.n_etats} -> ${d.n_etats_reduit} etats`;
    rendreSimplification();
  } catch (err) {
    $("#simp-etat").textContent = "";
    $("#simp-alertes").append(el("div", { class: "pb erreur" }, err.message));
  }
}

const SIMP_STATUTS = {
  retrouve: "retrouve",
  retire: "retire",
  perdu: "perdu",
  nouveau: "nouveau",
};

function rendreSimplification() {
  const d = E.res.simplification;
  if (!d) return;
  const al = $("#simp-alertes");
  al.textContent = "";
  const c = d.comptes;
  al.append(el("div", { class: "pb " + (d.sure ? "ok" : "avis") },
    `${d.sure ? "Aucune information de stabilite perdue." : "Information de stabilite perdue ou inventee :"} ` +
    `${d.n_etats} -> ${d.n_etats_reduit} etats ; max Re ${fmt(d.max_re, 4)} -> ${fmt(d.max_re_reduit, 4)} ; ` +
    `${c.retrouves} retrouve(s), ${c.retires} retire(s), ${c.perdus} perdu(s), ${c.nouveaux} nouveau(x) ; ` +
    `ecart max ${fmt(100 * d.ecart_max, 3)} % (${fmt(100 * d.ecart_max_em, 3)} % sous 5 Hz), dzeta max ${fmt(d.dzeta_max, 3)} pts`));
  for (const a of d.alertes) al.append(el("div", { class: "pb erreur" }, a));
  const zeta = (m) => (m.zeta ?? m.zeta_r);
  const lignesM = d.modes.filter((m) => (m.im ?? m.im_r) >= -1e-9).sort((a, b) => zeta(a) - zeta(b));
  tableau("#simp-modes",
    ["statut", "f [Hz]", "zeta [%]", "f simplifie [Hz]", "zeta simplifie [%]", "ecart [%]", "dzeta [pts]", "part retiree [%]", "elements"],
    lignesM.map((m) => ({
      m,
      cellules: [SIMP_STATUTS[m.statut], m.f ?? "", m.zeta ?? "", m.f_r ?? "", m.zeta_r ?? "",
        m.ecart == null ? "" : 100 * m.ecart, m.dzeta ?? "", 100 * m.part_retiree,
        m.elements.map((e) => `${e.nom} ${Math.round(100 * e.part)}%`).join(", ")],
    })),
    (l) => ((l.m.re ?? 0) > 0 || (l.m.re_r ?? 0) > 0 ? "instable" : l.m.critique ? "voisin" : l.m.statut === "retire" ? "doux" : ""));
  dessinerSimplification();
}

function dessinerSimplification() {
  const d = E.res.simplification;
  if (!d) return;
  const svg = $("#carte-simplification");
  svg.textContent = "";
  barreBornes("simp", dessinerSimplification);
  const xs = [], ys = [];
  for (const m of d.modes) {
    for (const [re, im] of [[m.re, m.im], [m.re_r, m.im_r]])
      if (re != null) { xs.push(re); ys.push(im, -im); }
  }
  const b = bornesPlan("simp", xs, ys, { zeroX: true, zeroY: true });
  const P = planLineaire(svg, b, { W: 1000, H: 480, ML: 70, MR: 18, MT: 30, MB: 44, titreX: "Re(lambda) [1/s]", titreY: "Im(lambda) [rad/s]" });
  isoAmortissement(P, false, [0.05, 0.1]);
  const dedans = (re, im) => re >= b.x0 && re <= b.x1 && im >= b.y0 && im <= b.y1;
  const croix = (x, y, col) => el("path", { d: `M${x - 5},${y - 5}l10,10M${x - 5},${y + 5}l10,-10`, stroke: col, "stroke-width": 2 });
  const bulle = (n, m) => infobulle(n, () =>
    `${m.statut}\n` + (m.re != null ? `complet : ${fmt(m.re, 4)} ${m.im >= 0 ? "+" : "-"} ${fmt(Math.abs(m.im), 4)}j, ${fmt(m.f, 3)} Hz, zeta ${fmt(m.zeta, 2)} %\n` : "") +
    (m.re_r != null ? `simplifie : ${fmt(m.re_r, 4)} ${m.im_r >= 0 ? "+" : "-"} ${fmt(Math.abs(m.im_r), 4)}j, ${fmt(m.f_r, 3)} Hz, zeta ${fmt(m.zeta_r, 2)} %\n` : "") +
    `${m.elements.map((e) => `${e.nom} ${Math.round(100 * e.part)} %`).join(", ")}`);
  for (const m of d.modes) {
    const conj = (im) => (Math.abs(im) > 1e-9 ? [1, -1] : [1]);
    if (m.statut === "retrouve") {
      for (const s of conj(m.im)) {
        if (dedans(m.re, s * m.im) || dedans(m.re_r, s * m.im_r))
          P.g.append(el("line", { x1: P.PX(m.re), y1: P.PY(s * m.im), x2: P.PX(m.re_r), y2: P.PY(s * m.im_r), stroke: "#1baf7a", "stroke-width": 1, opacity: 0.6 }));
        if (dedans(m.re_r, s * m.im_r)) {
          const x = P.PX(m.re_r), y = P.PY(s * m.im_r);
          const n = el("path", { d: `M${x},${y - 6}l6,6l-6,6l-6,-6z`, fill: "none", stroke: "#1baf7a", "stroke-width": 1.6 });
          bulle(n, m);
          P.g.append(n);
        }
      }
    }
    if (m.re != null) {
      for (const s of conj(m.im)) {
        if (!dedans(m.re, s * m.im)) continue;
        const x = P.PX(m.re), y = P.PY(s * m.im);
        const n = m.statut === "retrouve"
          ? el("circle", { cx: x, cy: y, r: 4, fill: "#1d2433" })
          : croix(x, y, m.statut === "perdu" ? "#d03b3b" : "#9aa1ad");
        bulle(n, m);
        P.g.append(n);
      }
    } else {
      for (const s of conj(m.im_r)) {
        if (!dedans(m.re_r, s * m.im_r)) continue;
        const n = el("rect", { x: P.PX(m.re_r) - 5, y: P.PY(s * m.im_r) - 5, width: 10, height: 10, fill: "none", stroke: "#7a4fa3", "stroke-width": 1.8 });
        bulle(n, m);
        P.g.append(n);
      }
    }
  }
  zoomPlan(svg, "simp", P);
  svg.dataset.nom = "simplification";
  svg.dataset.csv = "statut,Re,Im,f_Hz,zeta,Re_simplifie,Im_simplifie,f_simplifie_Hz,zeta_simplifie,part_retiree\n" +
    d.modes.map((m) => [m.statut, m.re, m.im, m.f, m.zeta, m.re_r, m.im_r, m.f_r, m.zeta_r, m.part_retiree].map((x) => x ?? "").join(",")).join("\n");
}

async function balayer() {
  const v = E.version;
  const a = parseFloat($("#bal-min").value), b = parseFloat($("#bal-max").value);
  const n = Math.max(2, Math.min(60, parseInt($("#bal-n").value, 10) || 20));
  const valeurs = Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
  $("#bal-etat").textContent = `calcul de ${n} points...`;
  $("#bal-critique").textContent = "";
  try {
    const d = await api("/api/balayage", corpsReseau({
      cible: cibleBalayage(),
      valeurs,
      f_max: null,
      n_modes: 100000, // TOUS les modes : l'utilisateur zoome lui-meme
    }));
    if (v !== E.version) return;
    E.res.balayage = d;
    E.lieux = null;
    const rates = d.pas.filter((p) => p.erreur).length;
    $("#bal-etat").textContent =
      `${n} points` + (rates ? `, dont ${rates} sans point de fonctionnement` : "");
    rendreBalayage();
  } catch (err) {
    $("#bal-etat").textContent = "";
    $("#carte-lieux").textContent = "";
    $("#bal-critique").textContent = "";
    $("#bal-critique").append(el("div", { class: "pb erreur" }, err.message));
  }
}

function rendreBalayage() {
  const d = E.res.balayage;
  if (!d) return;
  animerLieux(false);
  $("#bal-curseur-bloc").hidden = true;
  $("#bal-participations").textContent = "";
  const bons = d.pas.filter((p) => p.modes && p.modes.length);
  if (!bons.length) {
    $("#carte-lieux").textContent = "";
    $("#bal-critique").append(el("div", { class: "pb" }, "aucun point de fonctionnement sur cette plage"));
    return;
  }
  const nl = Math.min(...bons.map((p) => p.modes.length));
  E.lieux = { bons, nl };
  dessinerLieux();
  const r = $("#bal-curseur");
  r.max = bons.length - 1;
  r.value = 0;
  $("#bal-curseur-bloc").hidden = false;
  placerCurseurLieux(0);

  // --- marge de stabilite -------------------------------------------
  const g = $("#carte-marge");
  g.textContent = "";
  const W2 = 430, H2 = 300, ML2 = 58, MR2 = 14, MT2 = 14, MB2 = 38;
  const vx = bons.map((p) => p.valeur), vy = bons.map((p) => p.max_re);
  const ax = Math.min(...vx), bx = Math.max(...vx);
  const ay = Math.min(0, ...vy), by = Math.max(0, ...vy);
  const mx = (v2) => ML2 + ((v2 - ax) / (bx - ax || 1)) * (W2 - ML2 - MR2);
  const my = (v2) => H2 - MB2 - ((v2 - ay) / (by - ay || 1)) * (H2 - MT2 - MB2);
  g.append(el("line", { x1: ML2, y1: my(0), x2: W2 - MR2, y2: my(0), class: "axe" }));
  for (let k = 0; k <= 4; k++) {
    const v2 = ax + ((bx - ax) * k) / 4;
    g.append(el("text", { x: mx(v2), y: H2 - MB2 + 14, "text-anchor": "middle", class: "etiquette" }, fmt(v2, 4)));
  }
  for (let k = 0; k <= 4; k++) {
    const v2 = ay + ((by - ay) * k) / 4;
    g.append(el("line", { x1: ML2, y1: my(v2), x2: W2 - MR2, y2: my(v2), class: "grille" }));
    g.append(el("text", { x: ML2 - 6, y: my(v2) + 4, "text-anchor": "end", class: "etiquette" }, fmt(v2, 2)));
  }
  g.append(el("polyline", {
    points: bons.map((p) => `${mx(p.valeur)},${my(p.max_re)}`).join(" "),
    fill: "none", stroke: "#1f4e79", "stroke-width": 2,
  }));
  for (const p of bons) {
    const c = el("circle", { cx: mx(p.valeur), cy: my(p.max_re), r: 3, fill: p.max_re > 0 ? "#bf2617" : "#1f4e79" });
    infobulle(c, `${fmt(p.valeur, 5)} -> max Re = ${fmt(p.max_re, 4)}`);
    g.append(c);
  }
  g.append(el("text", { x: (ML2 + W2) / 2, y: H2 - 6, "text-anchor": "middle", class: "etiquette" }, cleCible()));
  g.append(el("text", { x: 10, y: MT2 - 2, class: "etiquette" }, "max Re [1/s]"));
  g.dataset.nom = "marge";
  g.dataset.csv = "valeur,max_re" + "\n" + bons.map((p) => [p.valeur, p.max_re].join(",")).join("\n");

  // --- la valeur critique, s'il y en a une ---------------------------
  const n = $("#bal-critique");
  n.textContent = "";
  let trouve = null;
  for (let i = 1; i < bons.length; i++) {
    const a0 = bons[i - 1], a1 = bons[i];
    if (a0.max_re < 0 && a1.max_re >= 0) {
      // interpolation lineaire du passage par zero : le pas du balayage est
      // bien plus grossier que la precision qu'on peut en tirer
      const x = a0.valeur + ((a1.valeur - a0.valeur) * (0 - a0.max_re)) / (a1.max_re - a0.max_re);
      const j = a1.modes.reduce((best, m, k) => (m.re > a1.modes[best].re ? k : best), 0);
      trouve = { x, f: a1.modes[j].f };
      break;
    }
  }
  if (trouve) {
    n.append(el("div", { class: "pb erreur" },
      `Valeur critique : ${cleCible()} = ${fmt(trouve.x, 5)} -- au-dela, le reseau est instable ` +
      `(mode a ${fmt(trouve.f, 2)} Hz).`));
  } else if (bons.every((p) => p.max_re < 0)) {
    n.append(el("div", { class: "pb" }, "stable sur toute la plage."));
  } else {
    n.append(el("div", { class: "pb erreur" }, "instable des le premier point de la plage."));
  }
}

// --- lieu des racines : plan complexe, echelle lineaire, tous les modes ---
// Un trait par mode (suivi d'un pas a l'autre), d'une seule couleur : avec
// des centaines de lieux, une couleur par lieu ne dirait plus rien. Le lieu
// qu'on survole dans le tableau se surligne.
function dessinerLieux() {
  const L = E.lieux;
  if (!L) return;
  const { bons, nl } = L;
  const svg = $("#carte-lieux");
  svg.textContent = "";
  barreBornes("bal", dessinerLieux);
  const tous = bons.flatMap((p) => p.modes.slice(0, nl));
  const b = bornesPlan("bal", tous.map((m) => m.re), tous.map((m) => m.im), { zeroX: true, zeroY: true });
  const P = planLineaire(svg, b, { W: 1000, H: 480, ML: 70, MR: 18, MT: 26, MB: 44, titreX: "Re(lambda) [1/s]", titreY: "Im(lambda) [rad/s]" });
  isoAmortissement(P, false, [0.05, 0.1]);
  L.traits = [];
  for (let j = 0; j < nl; j++) {
    const pts = bons.map((p) => `${P.PX(p.modes[j].re).toFixed(1)},${P.PY(p.modes[j].im).toFixed(1)}`);
    const t = el("polyline", { points: pts.join(" "), fill: "none", stroke: "#1f4e79", "stroke-width": 1.5, "stroke-opacity": 0.8 });
    P.g.append(t);
    L.traits.push(t);
    const d0 = bons[0].modes[j], dn = bons[bons.length - 1].modes[j];
    P.g.append(el("circle", { cx: P.PX(d0.re), cy: P.PY(d0.im), r: 3.5, fill: "#fff", stroke: "#1f4e79", "stroke-width": 1.5 }));
    const fin = el("circle", { cx: P.PX(dn.re), cy: P.PY(dn.im), r: 4, fill: dn.re > 0 ? "#bf2617" : "#1f4e79", stroke: "#fff", "stroke-width": 1 });
    infobulle(fin, () => `${fmt(dn.re, 4)} ${dn.im >= 0 ? "+" : "-"} ${fmt(Math.abs(dn.im), 4)}j` + "\n" +
      `${fmt(dn.f, 3)} Hz, zeta = ${fmt(dn.zeta, 2)} %` + "\n" + `a ${fmt(bons[bons.length - 1].valeur, 5)}`);
    P.g.append(fin);
  }
  L.P = P;
  L.curseur = el("g");
  P.g.append(L.curseur);
  zoomPlan(svg, "bal", P);
  svg.dataset.nom = "lieu_des_racines";
  svg.dataset.csv =
    "valeur,lieu,Re,Im,f_Hz,zeta_pct" + "\n" +
    bons.flatMap((p) => p.modes.slice(0, nl).map((m, j) => [p.valeur, j, m.re, m.im, m.f, m.zeta].join(","))).join("\n");
  placerCurseurLieux(+$("#bal-curseur").value || 0);
}

// --- parcourir le lieu ---------------------------------------------
// Le lieu montre OU passent les modes ; le curseur dit QUI oscille a chaque
// pas. Un mode de machine qui devient, a l'approche de l'instabilite, un
// mode de convertisseur ne se voit que la : le trait, lui, reste continu.
function placerCurseurLieux(i) {
  const L = E.lieux;
  if (!L || !L.P) return;
  i = Math.max(0, Math.min(L.bons.length - 1, i));
  const p = L.bons[i], P = L.P, b = P.b;
  // la position de TOUS les modes a ce pas, d'un seul chemin
  let chemin = "";
  for (const m of p.modes.slice(0, L.nl)) {
    const x = P.PX(m.re), y = P.PY(m.im);
    chemin += `M${(x + 5).toFixed(1)},${y.toFixed(1)}a5,5 0 1,0 -10,0a5,5 0 1,0 10,0`;
  }
  L.curseur.textContent = "";
  L.curseur.append(el("path", { d: chemin, fill: "#e08b2a", "fill-opacity": 0.25, stroke: "#e08b2a", "stroke-width": 1.5 }));

  // le tableau : les modes VISIBLES (dans les bornes), de frequence >= 0,
  // du moins amorti au plus amorti
  const tb = $("#bal-participations");
  tb.textContent = "";
  tb.append(el("thead", {}, el("tr", {},
    el("th", {}, "lieu"), el("th", {}, "Re"), el("th", {}, "Im"), el("th", {}, "f [Hz]"), el("th", {}, "zeta [%]"),
    el("th", { style: "text-align:left" }, "participations (par element)"),
    el("th", { style: "text-align:left" }, "etat dominant"))));
  const corps = el("tbody");
  const visibles = p.modes.slice(0, L.nl)
    .map((m, j) => ({ m, j }))
    .filter(({ m }) => m.im >= 0 && m.re >= b.x0 && m.re <= b.x1 && m.im >= b.y0 && m.im <= b.y1)
    .sort((u, v) => u.m.zeta - v.m.zeta)
    .slice(0, 20);
  for (const { m, j } of visibles) {
    const barres = el("td", { style: "text-align:left" });
    for (const [nom, part] of m.part || []) {
      barres.append(el("span", { class: "part" },
        el("span", { class: "part-barre", style: `width:${Math.round(part * 60)}px;background:#1f4e79` }),
        ` ${nom} ${Math.round(part * 100)} %`));
    }
    const tr = el("tr", {},
      el("td", { class: "mono" }, "#" + j),
      el("td", { class: "chiffre" }, fmt(m.re, 4)),
      el("td", { class: "chiffre" }, fmt(m.im, 4)),
      el("td", { class: "chiffre" }, fmt(m.f, 3)),
      el("td", { class: m.re > 0 ? "chiffre instable" : "chiffre" }, fmt(m.zeta, 2)),
      barres,
      el("td", { class: "mono", style: "text-align:left" }, m.etat || ""));
    const trait = L.traits[j];
    tr.addEventListener("mouseenter", () => {
      trait.setAttribute("stroke", "#e08b2a");
      trait.setAttribute("stroke-width", 3.5);
      trait.parentNode.append(trait); // au premier plan
    });
    tr.addEventListener("mouseleave", () => {
      trait.setAttribute("stroke", "#1f4e79");
      trait.setAttribute("stroke-width", 1.5);
    });
    corps.append(tr);
  }
  if (!visibles.length)
    corps.append(el("tr", {}, el("td", { colspan: 7, class: "doux" }, "aucun mode dans la fenetre : elargir les bornes (Tout voir)")));
  tb.append(corps);
  $("#bal-curseur").value = i;
  $("#bal-valeur").textContent = `${cleCible()} = ${fmt(p.valeur, 5)}` +
    (p.max_re > 0 ? "  -- instable" : "");
}

function animerLieux(marche) {
  clearInterval(E._animLieux);
  E._animLieux = null;
  $("#bal-anim").innerHTML = marche ? "&#10074;&#10074;" : "&#9654;";
  if (!marche || !E.lieux) return;
  const r = $("#bal-curseur");
  if (+r.value >= +r.max) placerCurseurLieux(0);
  E._animLieux = setInterval(() => {
    const i = +r.value + 1;
    if (!E.lieux || i > +r.max) return animerLieux(false);
    placerCurseurLieux(i);
  }, Math.max(80, 3000 / (+r.max + 1)));
}

function cleCible() {
  const e = elementBalayage(), q = parametreBalayage();
  return e && q ? `${e.nom} . ${q.nom}` : "";
}

// =====================================================================
// fichiers et export
// =====================================================================
/** Le reseau du constructeur, tel que le serveur l'attend : structure s'il a
 *  des zones, a plat sinon. */
function corpsFichier() {
  if (E.reseau) return corpsStructure();
  const tables = { ...E.tables, SIMPLIFICATIONS: E.simpAppliquer ? E.simp : {} };
  return { nom: E.tables.nom, tables, positions: E.positions };
}

// Telecharger un VRAI .g2egrid.json -- celui que le catalogue relit, et qu'on
// peut deposer dans g2egrid/reseaux/. Le format est l'affaire du serveur : le
// reecrire ici en ferait une seconde definition, qui deriverait.
async function enregistrer() {
  try {
    const d = await api("/api/fichiers/convertir", { reseau: corpsFichier() });
    const b = new Blob([JSON.stringify(d, null, 1)], { type: "application/json" });
    const a = el("a", { href: URL.createObjectURL(b), download: (d.nom || "reseau") + ".g2egrid.json" });
    document.body.append(a);
    a.click();
    a.remove();
    etat("telecharge : " + (d.nom || "reseau") + ".g2egrid.json");
  } catch (err) {
    etat("enregistrement impossible : " + err.message, true);
  }
}

async function enregistrerAuCatalogue() {
  try {
    const d = await api("/api/fichiers", { reseau: corpsFichier(), nom: $("#nom-reseau").value || null });
    etat(`enregistre dans le catalogue : ${d.nom} (${Math.round(d.octets / 1024)} ko)`);
    await rafraichirCatalogue(d.nom);
  } catch (err) {
    etat("enregistrement impossible : " + err.message, true);
  }
}

// Ouvrir n'importe quelle forme -- fichier natif, ancien JSON a plat -- :
// le serveur la ramene a la forme que le constructeur sait editer.
async function ouvrir(fichier) {
  try {
    const texte = await fichier.text();
    const ext = fichier.name.toLowerCase().split(".").pop();
    // un rapport d'etude porte ses tables : on le rouvre comme un reseau
    if (ext === "html" || ext === "htm") {
      const doc = new DOMParser().parseFromString(texte, "text/html");
      const n = doc.getElementById("g2egrid-tables");
      if (!n) throw new Error("ce fichier HTML n'est pas un rapport g2egrid");
      const d = JSON.parse(n.textContent);
      sortirZones();
      poser(d.tables, d.positions, d.nom);
      etat(`rapport rouvert : ${d.nom} (empreinte ${d.empreinte.slice(0, 12)}...)`);
      return;
    }
    // un fichier d'un autre outil : MATPOWER, PSS/E, ou le JSON de pandapower
    if (ext === "m" || ext === "raw" || (ext === "json" && texte.slice(0, 2000).includes('"pandapower'))) {
      const d = await api("/api/fichiers/externe", { nom: fichier.name, contenu: texte });
      sortirZones();
      poser(d.tables, d.positions, d.nom);
      const n = d.tables.Y_DER.length;
      etat(`importe : ${fichier.name} (${n} bus)` +
        (d.avertissements.length ? " -- " + d.avertissements.join(" ; ") : ""), d.avertissements.length > 0);
      return;
    }
    const brut = JSON.parse(texte);
    const d = await api("/api/fichiers/lire", { reseau: brut });
    if (d.instances && d.instances.length) entrerZones(d);
    else {
      sortirZones();
      poser(d.tables, d.positions, d.nom);
    }
    // Un fichier qui s'ouvre sans un seul bus n'est presque jamais voulu :
    // le dire, plutot que d'afficher un schema vide sous "ouvert".
    const n = (d.tables && d.tables.Y_DER ? d.tables.Y_DER.length : 0);
    if (!n) etat(`ouvert : ${fichier.name} -- mais il ne contient AUCUN bus`, true);
    else etat(`ouvert : ${fichier.name} (${n} bus)`);
  } catch (err) {
    etat("fichier illisible : " + err.message, true);
  }
}

// --- rapport d'etude ----------------------------------------------------
// Les figures deja calculees partent avec la demande, en SVG : le rapport
// montre ce qu'on a vu. Le reste (load flow, spectre, modes critiques) est
// recalcule par le serveur, sur les tables dont le rapport donne l'empreinte.
function figuresCalculees() {
  const ser = new XMLSerializer();
  const f = {};
  const prendre = (legende, sel, condition = true) => {
    const n = $(sel);
    if (!condition || !n || !n.childElementCount) return;
    const c = n.cloneNode(true);
    c.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    c.removeAttribute("id");
    f[legende] = ser.serializeToString(c);
  };
  prendre("Schema du reseau" + (E.vue && E.vue !== "transport" ? ` (zone ${E.vue})` : ""), "#canevas");
  prendre("Point de fonctionnement", "#canevas-lf", !!E.res.lf);
  prendre("Carte des modes", "#carte-modes", !!E.res.modal);
  prendre("Mode sur le schema", "#schema-mode", !!E.res.forme);
  prendre("Continuation : max Re", "#carte-cont", !!E.res.continuation);
  prendre("Continuation : lieu des modes", "#carte-lieu", !!E.res.continuation);
  prendre("Boucle mineure : lieux caracteristiques", "#carte-nyquist", !!E.res.impedance);
  prendre("Boucle mineure : impedances", "#carte-zbode", !!E.res.impedance);
  return f;
}

async function rapportEtude() {
  try {
    dessiner();
    const d = await lancerTache("rapport", corpsReseau({ figures: figuresCalculees() }));
    telecharger(d.html, "text/html", d.fichier);
    etat(`${d.fichier} telecharge (empreinte ${d.empreinte.slice(0, 12)}...)`);
  } catch (err) {
    etat("rapport impossible : " + err.message, true);
  }
}

async function exporter() {
  try {
    const d = await api("/api/exporter", corpsReseau());
    $("#sortie-export").value = d.python;
    $("#export-etat").textContent = "cellule de donnees : collable telle quelle dans un notebook";
  } catch (err) {
    $("#sortie-export").value = "ECHEC : " + err.message;
  }
}

// Les exports "lourds" : le texte s'affiche ET se telecharge, et ce que
// l'export perd -- shunts en commentaire dans le .m... -- est DIT.
async function exporterVers(quoi) {
  try {
    const d = await api("/api/exporter/" + quoi, corpsReseau());
    $("#sortie-export").value = d.contenu;
    telecharger(d.contenu, d.type || "text/plain", d.fichier);
    $("#export-etat").textContent = d.avertissements.length
      ? "ATTENTION : " + d.avertissements.join(" ; ")
      : `${d.fichier} telecharge`;
    $("#export-etat").style.color = d.avertissements.length ? "var(--erreur)" : "var(--doux)";
  } catch (err) {
    $("#export-etat").textContent = "impossible : " + err.message;
    $("#export-etat").style.color = "var(--erreur)";
  }
}

async function produireCode() {
  try {
    const nom = $("#code-catalogue").checked ? ($("#nom-reseau").value || null) : null;
    const d = await api("/api/code", { reseau: corpsFichier(), nom });
    $("#code-python").value = d.python;
  } catch (err) {
    $("#code-python").value = "ECHEC : " + err.message;
  }
}

// =====================================================================
// onglets et demarrage
// =====================================================================
// =====================================================================
// documentation
// =====================================================================
// La doc Sphinx est servie par l'atelier sous /doc ; si elle n'est pas encore
// construite, l'onglet propose de la construire (Sphinx doit etre installe).
async function etatDoc() {
  const d = await api("/api/doc");
  const n = $("#doc-etat");
  if (d.construction === "en cours") n.textContent = "construction en cours (une minute environ)...";
  else if (d.construction === "echec") n.textContent = "la construction a echoue : " + d.message.slice(-300);
  else n.textContent = d.construite ? "" : "documentation pas encore construite";
  return d;
}

async function ouvrirDoc(page = "index.html") {
  // marque AVANT d'activer l'onglet : activerOnglet("doc") appelle ouvrirDoc()
  // tant que rien n'a ete demande, et le cadre ne recoit son adresse qu'apres
  // l'attente du serveur -- sans cette marque, les deux fonctions s'appelaient
  // l'une l'autre jusqu'a epuiser la pile, au premier clic sur l'onglet
  E.docDemandee = true;
  activerOnglet("doc");
  const d = await etatDoc();
  if (d.construite) {
    $("#cadre-doc").src = "/doc/" + page;
    $("#lien-doc").href = "/doc/" + page;
  } else {
    $("#cadre-doc").srcdoc =
      "<p style='font-family:sans-serif;padding:20px'>La documentation n'est pas encore construite. " +
      "Bouton <b>Construire</b> ci-dessus, ou en ligne de commande : " +
      "<code>uv run --group docs sphinx-build -b html docs docs/_build/html</code></p>";
  }
}

async function construireDoc() {
  try {
    await api("/api/doc/construire", {});
  } catch (err) {
    $("#doc-etat").textContent = err.message;
    return;
  }
  const suivre = async () => {
    const d = await etatDoc();
    if (d.construction === "en cours") setTimeout(suivre, 3000);
    else if (d.construite) ouvrirDoc();
  };
  suivre();
}

// =====================================================================
// onglet Limites et impedances
// =====================================================================
const TEINTES = ["#1f4e79", "#bf2617", "#2a7a4b", "#e08b2a", "#6b3fa0", "#0f7d8c", "#8a5a2b", "#c2185b"];

/** Un trace en lignes : series = [{xs, ys, couleur, largeur, pointille, nom}],
 *  verticales = [{x, couleur, texte}], marqueurs = [{x, y, couleur, texte}]. */
function tracerCourbes(svg, prefixe, series, o) {
  const W = 1000, H = +svg.viewBox.baseVal.height || 380;
  const tx = o.logX ? (x) => Math.log10(x) : (x) => x;
  const xs = series.flatMap((s) => s.xs.map(tx)).filter(Number.isFinite);
  const ys = series.flatMap((s) => s.ys).filter(Number.isFinite);
  for (const m of o.marqueurs || []) {
    xs.push(tx(m.x));
    ys.push(m.y);
  }
  const rendre = () => {
    svg.textContent = "";
    const b = bornesPlan(prefixe, xs, ys, { zeroY: !!o.zeroY, zeroX: !!o.zeroX });
    const P = planLineaire(svg, b, { W, H, ML: 70, MR: 20, MT: 26, MB: 42, titreX: o.titreX, titreY: o.titreY });
    for (const [k, s] of series.entries()) {
      const pts = s.xs.map((x, i) => [tx(x), s.ys[i]]).filter(([x, y]) => Number.isFinite(x) && Number.isFinite(y));
      if (!pts.length) continue;
      const couleur = s.couleur || TEINTES[k % TEINTES.length];
      const ligne = el("polyline", {
        points: pts.map(([x, y]) => `${P.PX(x)},${P.PY(y)}`).join(" "),
        fill: "none", stroke: couleur, "stroke-width": s.largeur || 1.8,
        "stroke-dasharray": s.pointille ? "5 4" : null,
      });
      if (s.nom) infobulle(ligne, s.nom);
      P.g.append(ligne);
    }
    for (const v of o.verticales || []) {
      const X = P.PX(tx(v.x));
      P.g.append(el("line", { x1: X, y1: 0, x2: X, y2: H, stroke: v.couleur || "#bf2617", "stroke-width": 1.5, "stroke-dasharray": "6 3" }));
      const droite = X > W - 140; // pres du bord : le texte passe a gauche du trait
      if (v.texte) P.g.append(el("text", { x: droite ? X - 4 : X + 4, y: 40, "text-anchor": droite ? "end" : "start", class: "etiquette", fill: v.couleur || "#bf2617" }, v.texte));
    }
    for (const m of o.marqueurs || []) {
      const X = P.PX(tx(m.x)), Y = P.PY(m.y);
      const c = el("polygon", { points: `${X},${Y - 7} ${X + 7},${Y} ${X},${Y + 7} ${X - 7},${Y}`, fill: m.couleur || "#bf2617", stroke: "#fff" });
      if (m.texte) infobulle(c, m.texte);
      P.g.append(c);
    }
    if (o.apres) o.apres(P);
    zoomPlan(svg, prefixe, P);
  };
  barreBornes(prefixe, rendre);
  rendre();
}

function listeBifurcations(n, bifs) {
  n.textContent = "";
  for (const b of bifs)
    n.append(el("div", { class: "pb " + (b.sens > 0 ? "erreur" : "avis") }, b.texte));
  if (!bifs.length)
    n.append(el("div", { class: "pb ok" }, "aucune bifurcation sur l'intervalle"));
}

async function chargerElements() {
  if (E._elementsVersion === E.version && E._elements) return E._elements;
  const d = await api("/api/elements", corpsReseau());
  E._elements = d;
  E._elementsVersion = E.version;
  const remplir = (sel, liste) => {
    const garde = sel.value;
    sel.textContent = "";
    for (const x of liste) sel.append(el("option", { value: x.nom }, `${x.nom} (${x.classe})`));
    if (liste.some((x) => x.nom === garde)) sel.value = garde;
  };
  remplir($("#cont-element"), d);
  remplir($("#imp-element"), d.filter((x) => x.convertisseur));
  remplirParametresContinuation();
  return d;
}

// (nom distinct de remplirParametres, celle de l'onglet Variations : deux
// declarations de meme nom, la seconde ecrase la premiere et la liste des
// parametres a balayer restait vide)
function remplirParametresContinuation() {
  const x = (E._elements || []).find((e) => e.nom === $("#cont-element").value);
  const sel = $("#cont-param");
  const garde = sel.value;
  sel.textContent = "";
  if (!x) return;
  for (const [k, v] of Object.entries(x.parametres)) sel.append(el("option", { value: k }, `${k} = ${fmt(v, 4)}`));
  if (garde in x.parametres) sel.value = garde;
  bornesParametre();
}

function bornesParametre() {
  const x = (E._elements || []).find((e) => e.nom === $("#cont-element").value);
  if (!x) return;
  const v = x.parametres[$("#cont-param").value];
  if (v === undefined) return;
  const base = Math.abs(v) || 1;
  $("#cont-de").value = +(base / 3).toPrecision(4);
  $("#cont-a").value = +(base * 3).toPrecision(4);
}

/** Les stations en statisme Vdc-P du reseau : bus alternatif -> kdroop. */
function stationsStatisme() {
  const out = {};
  for (const r of lignes("T_VSCs"))
    if (String(r[3]) === "droop") out[+r[1]] = +(((E.tables.PARAMS_PAR_BUS || {})[r[1]] || {}).kdroop || 0);
  return out;
}

function remplirStatismes() {
  const sel = $("#stat-stations");
  const avant = sel.value;
  sel.textContent = "";
  const st = stationsStatisme();
  const bus = Object.keys(st);
  sel.append(el("option", { value: "" }, bus.length ? `toutes (${bus.length})` : "aucune station en statisme"));
  for (const b of bus) sel.append(el("option", { value: b }, `bus ${b} (k = ${fmt(100 * st[b], 3)} %)`));
  if ([...sel.options].some((o) => o.value === avant)) sel.value = avant;
}

async function balayerStatismes() {
  const v = E.version;
  $("#stat-etat").textContent = "calcul...";
  $("#stat-bifs").textContent = "";
  try {
    const choix = $("#stat-stations").value;
    const d = await lancerTache("statismes", corpsReseau({
      debut: parseFloat($("#stat-de").value) || 0.25,
      fin: parseFloat($("#stat-a").value) || 4,
      pas: parseInt($("#stat-pas").value) || 9,
      stations: choix ? [+choix] : null,
    }));
    if (v !== E.version) return;
    E.res.statismes = d;
    $("#stat-etat").textContent = `stations ${d.stations.join(", ")} ; ${d.resume.split("\n")[0]}`;
    const verticales = [{ x: 1, texte: "statismes du reseau" }];
    const enDroop = new Set(d.stations.map(String));
    const sp = Object.entries(d.p_mw).map(([b, ys], k) => ({
      xs: d.mu, ys, nom: `bus ${b}${enDroop.has(b) ? " (statisme)" : ""}`,
      couleur: TEINTES[k % TEINTES.length], largeur: enDroop.has(b) ? 2.4 : 1.2, pointille: !enDroop.has(b),
    }));
    tracerCourbes($("#carte-statp"), "statp", sp, { titreX: "facteur sur k (log10)", titreY: "P [MW]", logX: true, verticales });
    const sv = Object.entries(d.vdc).map(([b, ys], k) => ({ xs: d.mu, ys, nom: `noeud ${b}`, couleur: TEINTES[k % TEINTES.length], largeur: 1.4 }));
    tracerCourbes($("#carte-statv"), "statv", sv, { titreX: "facteur sur k (log10)", titreY: "Vdc [pu]", logX: true, verticales });
    tracerCourbes($("#carte-statre"), "statre", [{ xs: d.mu, ys: d.max_re, nom: "max Re(lambda)", largeur: 2.2 }], {
      titreX: "facteur sur k (log10)", titreY: "max Re [1/s]", logX: true, zeroY: true,
      verticales: [...verticales, ...d.bifurcations.map((b) => ({ x: b.mu, texte: `${b.genre} ${fmt(b.mu, 4)}` }))],
    });
    listeBifurcations($("#stat-bifs"), d.bifurcations);
  } catch (err) {
    $("#stat-etat").textContent = "";
    $("#stat-bifs").append(el("div", { class: "pb erreur" }, err.message));
  }
}

async function suivreParametre() {
  const v = E.version;
  $("#cont-etat").textContent = "calcul...";
  try {
    const d = await lancerTache("continuation", corpsReseau({
      element: $("#cont-element").value,
      parametre: $("#cont-param").value,
      debut: parseFloat($("#cont-de").value),
      fin: parseFloat($("#cont-a").value),
      pas: parseInt($("#cont-pas").value) || 15,
      log: $("#cont-log").checked,
      equilibre: $("#cont-eq").checked,
    }));
    if (v !== E.version) return;
    E.res.continuation = d;
    const nom = `${$("#cont-element").value}.${$("#cont-param").value}`;
    $("#cont-etat").textContent = d.resume.split("\n")[0];
    const log = $("#cont-log").checked;
    tracerCourbes($("#carte-cont"), "cont", [{ xs: d.mu, ys: d.max_re, nom: "max Re(lambda)", largeur: 2.2 }], {
      titreX: nom + (log ? " (log10)" : ""), titreY: "max Re(lambda) [1/s]", logX: log, zeroY: true,
      verticales: d.bifurcations.map((b) => ({ x: b.mu, texte: `${b.genre} ${fmt(b.mu, 4)}` })),
    });
    // le lieu des modes : du clair (debut) au fonce (fin)
    const series = [];
    for (const [k, t] of d.trajectoires.entries()) {
      const base = TEINTES[k % TEINTES.length];
      series.push({ xs: t.re, ys: t.im.map((y) => y / (2 * Math.PI)), couleur: base, largeur: 1.6, nom: `mode ${k + 1}` });
    }
    tracerCourbes($("#carte-lieu"), "lieu", series, {
      titreX: "Re(lambda) [1/s]", titreY: "f [Hz]", zeroX: true,
      apres: (P) => {
        for (const [k, t] of d.trajectoires.entries()) {
          const n = t.re.length;
          t.re.forEach((x, i) => {
            const c = el("circle", { cx: P.PX(x), cy: P.PY(t.im[i] / (2 * Math.PI)), r: 3, fill: TEINTES[k % TEINTES.length], "fill-opacity": 0.25 + (0.75 * i) / Math.max(1, n - 1) });
            infobulle(c, `${nom} = ${fmt(d.mu[i], 4)}` + "\n" + `lambda = ${fmt(x, 4)} + ${fmt(t.im[i], 4)}j`);
            P.g.append(c);
          });
        }
      },
    });
    listeBifurcations($("#cont-bifs"), d.bifurcations);
  } catch (err) {
    $("#cont-etat").textContent = "";
    $("#cont-bifs").textContent = "";
    $("#cont-bifs").append(el("div", { class: "pb erreur" }, err.message));
  }
}

/** Les lieux caracteristiques, sur tout l'axe (f < 0 : les conjugues). En
 *  echelle compressee, un point de rayon r est trace a log10(1 + r) : le
 *  point critique -1 tombe a -0,301, et les excursions lointaines des
 *  resonances ne l'ecrasent plus. Le nombre de tours n'en change pas. */
function rendreNyquist() {
  const d = E.res.impedance;
  if (!d) return;
  const comp = $("#imp-compresse").checked;
  const c = (x, y) => {
    if (!comp) return [x, y];
    const r = Math.hypot(x, y);
    const k = r > 0 ? Math.log10(1 + r) / r : 1;
    return [k * x, k * y];
  };
  const series = d.lieux.flatMap((l, k) => {
    const pts = l.re.map((x, i) => c(x, l.im[i]));
    return [
      { xs: pts.map((p) => p[0]), ys: pts.map((p) => p[1]), couleur: TEINTES[k], nom: `lieu ${k + 1}, f > 0` },
      { xs: pts.map((p) => p[0]), ys: pts.map((p) => -p[1]), couleur: TEINTES[k], pointille: true, nom: `lieu ${k + 1}, f < 0` },
    ];
  });
  const [mx, my] = c(-1, 0);
  tracerCourbes($("#carte-nyquist"), "nyq", series, {
    titreX: comp ? "Re (rayon en log10(1 + r))" : "Re", titreY: "Im", zeroX: true, zeroY: true,
    marqueurs: [{ x: mx, y: my, couleur: "#1c2330", texte: "point critique -1" }],
  });
}

async function boucleMineure() {
  const v = E.version;
  $("#imp-etat").textContent = "calcul...";
  try {
    const d = await lancerTache("impedance", corpsReseau({
      element: $("#imp-element").value,
      f_min: parseFloat($("#imp-fmin").value) || 0.1,
      f_max: parseFloat($("#imp-fmax").value) || 5000,
    }));
    if (v !== E.version) return;
    E.res.impedance = d;
    $("#imp-etat").textContent = `${d.element} au bus ${d.bus}, forme ${d.forme}`;
    const n = $("#imp-verdict");
    n.textContent = "";
    const genre = !d.applicable ? "avis" : d.stable ? "ok" : "erreur";
    const verdict = !d.applicable
      ? `boucle ouverte instable (P = ${d.instables + d.tours}) : le nombre de tours ne se lit pas seul`
      : d.stable ? "critere satisfait : aucun tour autour de -1" : `${-d.tours} tour(s) autour de -1 : ${d.instables} mode(s) instable(s)`;
    n.append(el("div", { class: "pb " + genre }, `${verdict} ; marge |1 + lambda(L)| = ${fmt(d.marge, 3)} a ${fmt(d.f_marge, 1)} Hz`));
    rendreNyquist();
    const dB = (z) => 20 * Math.log10(Math.max(z, 1e-12));
    tracerCourbes($("#carte-zbode"), "zb", [
      { xs: d.f, ys: d.z_reste.map(dB), couleur: "#1f4e79", nom: "|Z| du reste du reseau", largeur: 2 },
      { xs: d.f, ys: d.z_equipement.map(dB), couleur: "#2a7a4b", nom: `|Z| de ${d.element}`, largeur: 2 },
    ], {
      titreX: "f [Hz] (log10)", titreY: "|Z| [dB pu]", logX: true,
      apres: (P) => {
        for (const [a, b] of d.non_passif) {
          const x0 = P.PX(Math.log10(a)), x1 = P.PX(Math.log10(b));
          P.g.append(el("rect", { x: Math.min(x0, x1), y: 0, width: Math.max(2, Math.abs(x1 - x0)), height: 2000, fill: "#bf2617", "fill-opacity": 0.08 }));
        }
      },
    });
  } catch (err) {
    $("#imp-etat").textContent = "";
    $("#imp-verdict").textContent = "";
    $("#imp-verdict").append(el("div", { class: "pb erreur" }, err.message));
  }
}

// =====================================================================
// evenements : les grandeurs electriques apres un evenement
// =====================================================================
/* « Que devient le reseau apres ceci ? » Un scenario d'evenements
 * transitoires (ouverture, blocage, defaut, charge, saut de phase) ou de
 * consigne (echelon, rampe), des relais de surtension sur les
 * convertisseurs, et, pour CHAQUE bus et CHAQUE unite : l'avant, les
 * extremes du transitoire, le regime atteint. La cascade quasi statique pose
 * la meme question en regimes permanents successifs, sans dynamique. */
const EVT = { liste: [], res: null, cascade: null, bus: new Set(), unites: new Set() };
const EVT_CONSIGNES = ["echelon", "rampe"];
//: les valeurs que porte chaque genre : [cle, etiquette, valeur par defaut]
const EVT_CHAMPS = {
  echelon: [["amplitude", "Amplitude [pu]", 0.05]],
  rampe: [["amplitude", "Amplitude [pu]", 0.05], ["duree", "Duree [s]", 1]],
  charge: [["dp", "dP (relatif)", 0.1], ["dq", "dQ (relatif)", 0]],
  saut_de_phase: [["angle", "Angle [deg]", 10]],
  defaut: [["R", "R [pu]", 1e-3], ["X", "X [pu]", 1e-2]],
};

async function chargerEvenements() {
  if (!E.res.signaux) await chargerSignaux();
  remplirCiblesEvt(true);
}

function remplirCiblesEvt(nouveauGenre = false) {
  const genre = $("#evt-genre").value;
  const consigne = EVT_CONSIGNES.includes(genre);
  const toutes = consigne ? (E.res.signaux && E.res.signaux.entrees) || [] : ciblesManoeuvre(genre);
  const filtre = $("#evt-filtre").value.trim().toLowerCase();
  const s = $("#evt-cibles");
  const garde = new Set(Array.from(s.selectedOptions).map((o) => o.value));
  s.textContent = "";
  for (const c of toutes) {
    if (filtre && !c.toLowerCase().includes(filtre)) continue;
    const o = el("option", { value: c }, c);
    o.selected = garde.has(c);
    s.append(o);
  }
  $("#evt-cibles-lbl").textContent = consigne ? "Consigne(s)" : "Element(s)";
  const champs = EVT_CHAMPS[genre] || [];
  [["#evt-champ-valeur", "#evt-valeur-lbl", "#evt-valeur"], ["#evt-champ-valeur2", "#evt-valeur2-lbl", "#evt-valeur2"]].forEach(
    ([bloc, lbl, champ], k) => {
      $(bloc).hidden = !champs[k];
      if (!champs[k]) return;
      $(lbl).textContent = champs[k][1];
      if (nouveauGenre) $(champ).value = champs[k][2];
    },
  );
}

function ajouterEvt() {
  const genre = $("#evt-genre").value;
  const instant = Math.max(0, parseFloat($("#evt-instant").value) || 0);
  const cibles = Array.from($("#evt-cibles").selectedOptions).map((o) => o.value);
  if (!cibles.length) {
    $("#evt-etat").textContent = "choisir au moins un element (ou une consigne)";
    return;
  }
  const valeurs = {};
  (EVT_CHAMPS[genre] || []).forEach(([cle], k) => (valeurs[cle] = parseFloat($(k ? "#evt-valeur2" : "#evt-valeur").value) || 0));
  for (const c of cibles) {
    const ev = { instant, genre, ...valeurs };
    if (EVT_CONSIGNES.includes(genre)) ev.entree = c;
    else ev.cible = c;
    EVT.liste.push(ev);
  }
  rendreListeEvt();
}

function texteValeurEvt(ev) {
  return (EVT_CHAMPS[ev.genre] || []).map(([cle, lbl]) => `${lbl} ${fmt(ev[cle], 4)}`).join(" ; ");
}

function rendreListeEvt() {
  EVT.liste.sort((a, b) => a.instant - b.instant);
  const tb = $("#evt-table tbody");
  tb.textContent = "";
  EVT.liste.forEach((ev, i) => {
    const t = el("input", { type: "number", value: ev.instant, step: 0.1, min: 0, style: "width:70px" });
    t.addEventListener("change", () => {
      ev.instant = Math.max(0, parseFloat(t.value) || 0);
      rendreListeEvt();
    });
    const sup = el("button", { title: "Retirer cet evenement" }, "x");
    sup.addEventListener("click", () => {
      EVT.liste.splice(i, 1);
      rendreListeEvt();
    });
    tb.append(
      el("tr", {},
        el("td", {}, t),
        el("td", {}, ev.genre.replaceAll("_", " ")),
        el("td", { class: "mono" }, ev.cible || ev.entree),
        el("td", { class: "doux" }, texteValeurEvt(ev)),
        el("td", {}, sup)),
    );
  });
}

function relaisEvt() {
  if (!$("#evt-relais").checked) return null;
  const familles = [];
  if ($("#evt-fam-gfl").checked) familles.push("GFL");
  if ($("#evt-fam-gfm").checked) familles.push("GFM");
  const r = $("#evt-retombee").value.trim();
  return {
    v_max: parseFloat($("#evt-vmax").value) || 1.1,
    retard: Math.max(0, parseFloat($("#evt-retard").value) || 0),
    retombee: r === "" ? null : parseFloat(r),
    action: $("#evt-action").value,
    familles: familles.length ? familles : ["GFL"],
  };
}

async function simulerEvt() {
  const etat = $("#evt-etat");
  if (!EVT.liste.length) {
    etat.textContent = "le scenario est vide : ajouter au moins un evenement";
    return;
  }
  const v = E.version;
  const bouton = $("#btn-evt-simuler");
  bouton.disabled = true;
  etat.textContent = "simulation...";
  const t0 = performance.now();
  try {
    const pas = parseFloat($("#evt-pas").value);
    const d = await lancerTache("apres_evenement", corpsReseau({
      evenements: EVT.liste,
      duree: parseFloat($("#evt-duree").value) || 5,
      n: parseInt($("#evt-n").value, 10) || 501,
      pas: Number.isFinite(pas) && pas > 0 ? pas / 1000 : null,
      limiteurs: $("#evt-limiteurs").checked,
      relais: relaisEvt(),
      protections_dc: $("#evt-prot-dc").checked,
      v_bas: parseFloat($("#evt-vbas").value) || 0.9,
      v_haut: parseFloat($("#evt-vhaut").value) || 1.1,
    }));
    if (v !== E.version) return;
    EVT.res = d;
    choisirTracesEvt(d);
    etat.textContent = `simule en ${((performance.now() - t0) / 1000).toFixed(1)} s`;
    rendreEvt();
  } catch (e) {
    etat.textContent = e.message;
  } finally {
    bouton.disabled = false;
  }
}

// Les courbes tracees d'emblee : les bus qui ont le plus bouge et le plus
// haut, les unites sorties et celles dont le reactif a le plus change.
function choisirTracesEvt(d) {
  const bus = d.bus.filter((b) => b.V0 !== null);
  const ecart = [...bus].sort((a, b) => Math.abs(b.dV ?? 0) - Math.abs(a.dV ?? 0)).slice(0, 5);
  const haut = [...bus].sort((a, b) => b.Vmax - a.Vmax)[0];
  EVT.bus = new Set([...ecart, haut].filter(Boolean).map((b) => b.bus));
  const u = [...d.unites].sort((a, b) => Math.abs(b.dQ ?? 0) - Math.abs(a.dQ ?? 0));
  EVT.unites = new Set(u.slice(0, 6).map((x) => x.nom));
}

function rendreEvt() {
  const d = EVT.res;
  if (!d) return;
  const r = d.resume;
  const tuile = (val, txt) => el("div", {}, el("b", {}, val), el("span", {}, txt));
  const n = $("#evt-resume");
  n.className = "resume-evt";
  n.textContent = "";
  n.append(
    tuile(r.v_max ? `${fmt(r.v_max.V, 4)} pu` : "--", r.v_max ? `V max : bus ${r.v_max.bus} a ${fmt(r.v_max.t, 3)} s` : "V max"),
    tuile(r.v_min ? `${fmt(r.v_min.V, 4)} pu` : "--", r.v_min ? `V min : bus ${r.v_min.bus} a ${fmt(r.v_min.t, 3)} s` : "V min"),
    tuile(String(r.n_hors_bornes), `bus sortis de [${d.bornes[0]}, ${d.bornes[1]}] pu`),
    tuile(String(r.n_declenchements), "declenchements de relais"),
    tuile(String(r.n_sorties), "unites sorties en tout"),
    tuile(d.frequence.min !== null ? `${fmt(d.frequence.min, 3)} Hz` : "--", "frequence min"),
    tuile(d.regime_lf.ok ? "calcule" : "--", d.regime_lf.ok ? "regime du load flow apres les ouvertures" : `load flow : ${d.regime_lf.raison}`),
  );
  $("#evt-rangee-v").hidden = false;
  $("#evt-rangee-u").hidden = false;
  $("#evt-rangee-f").hidden = false;
  rendreTensionsEvt();
  rendreUnitesEvt();
  rendreFrequenceEvt();
  rendreRelaisEvt();
}

// Traits verticaux : les instants du scenario (gris), les declenchements (rouge).
function verticalesEvt(d) {
  const out = [];
  for (const t of new Set(d.evenements.map((e) => e.instant))) out.push({ x: t, couleur: "#8a93a5" });
  for (const x of d.declenchements) out.push({ x: x.t, couleur: "#bf2617", texte: x.protection.split(":").pop() });
  return out;
}

// Un trait horizontal en tirets, qu'on ne trace que s'il est pres des
// courbes : une borne a 0,9 pu ecraserait une excursion de 5 mpu.
function horizontalesEvt(d, series, lignes) {
  const ys = series.flatMap((s) => s.ys).filter(Number.isFinite);
  if (!ys.length) return;
  const [a, b] = minMax(ys);
  const marge = Math.max(0.02, 0.5 * (b - a));
  for (const [y, couleur, nom] of lignes) {
    if (y < a - marge || y > b + marge) continue;
    series.push({ xs: [d.t[0], d.t[d.t.length - 1]], ys: [y, y], couleur, pointille: true, largeur: 1.2, nom });
  }
}

const pastilleEvt = (k) =>
  k >= 0 ? el("span", { style: `display:inline-block;width:10px;height:10px;border-radius:2px;background:${TEINTES[k % TEINTES.length]}` }) : "";

function tableEvt(cible, entetes, lignes_) {
  const n = $(cible);
  n.textContent = "";
  const t = el("table");
  t.append(el("thead", {}, el("tr", {}, ...entetes.map((x) => el("th", {}, x)))));
  const tb = el("tbody");
  for (const l of lignes_) tb.append(l);
  t.append(tb);
  n.append(t);
}

function celluleEvt(v, n = 4, classe = "") {
  return el("td", { class: "chiffre " + classe }, v === null || v === undefined ? "--" : fmt(v, n));
}

function rendreTensionsEvt() {
  const d = EVT.res;
  const [vb, vh] = d.bornes;
  const choisis = [...EVT.bus];
  const series = choisis.map((b, k) => ({
    xs: d.t, ys: d.series.V[b].map((v) => v ?? NaN), nom: `bus ${b}`, couleur: TEINTES[k % TEINTES.length],
  }));
  const lignes_ = [[vh, "#8a93a5", `borne haute ${vh} pu`], [vb, "#8a93a5", `borne basse ${vb} pu`]];
  if (d.relais && Math.abs(d.relais.v_max - vh) > 1e-9) lignes_.push([d.relais.v_max, "#bf2617", `seuil des relais ${d.relais.v_max} pu`]);
  horizontalesEvt(d, series, lignes_);
  tracerCourbes($("#evt-courbes-v"), "evtv", series, { titreX: "temps [s]", titreY: "|V| [pu]", verticales: verticalesEvt(d) });
  $("#evt-v-titre").textContent = `${choisis.length} bus traces sur ${d.bus.length}`;

  const tri = $("#evt-tri").value;
  const filtre = $("#evt-filtre-bus").value.trim();
  let l = d.bus.filter((b) => b.V0 !== null);
  if (filtre) l = l.filter((b) => String(b.bus).startsWith(filtre));
  if ($("#evt-hors-seuls").checked) l = l.filter((b) => b.t_hors > 0);
  const cle = {
    dV: (b) => -Math.abs(b.dV ?? 0), Vmax: (b) => -b.Vmax, Vmin: (b) => b.Vmin, t_hors: (b) => -b.t_hors,
    ecart_lf: (b) => -Math.abs(b.ecart_lf ?? 0), bus: (b) => b.bus,
  }[tri];
  l.sort((a, b) => cle(a) - cle(b));
  const hors = (v) => (v !== null && (v > vh || v < vb) ? "v-hors" : "");
  const mpu = (v) => (v === null || v === undefined ? null : 1000 * v);
  tableEvt("#evt-table-bus",
    ["", "bus", "Un [kV]", "V0", "V min", "a [s]", "V max", "a [s]", "V fin", "dV [mpu]", "hors bornes [s]", "V lf", "fin - lf [mpu]"],
    l.map((b) => {
      const k = choisis.indexOf(b.bus);
      const tr = el("tr", { class: [k >= 0 ? "choisi" : "", b.t_hors > 0 ? "hors-bornes" : ""].join(" ") },
        el("td", {}, pastilleEvt(k)), el("td", { class: "mono" }, String(b.bus)), celluleEvt(b.Un_kV, 1),
        celluleEvt(b.V0), celluleEvt(b.Vmin, 4, hors(b.Vmin)), celluleEvt(b.t_min, 3),
        celluleEvt(b.Vmax, 4, hors(b.Vmax)), celluleEvt(b.t_max, 3), celluleEvt(b.Vfin, 4, hors(b.Vfin)),
        celluleEvt(mpu(b.dV), 2), celluleEvt(b.t_hors, 3), celluleEvt(b.V_lf), celluleEvt(mpu(b.ecart_lf), 2));
      tr.addEventListener("click", () => {
        if (EVT.bus.has(b.bus)) EVT.bus.delete(b.bus);
        else EVT.bus.add(b.bus);
        rendreTensionsEvt();
      });
      return tr;
    }));
}

function rendreUnitesEvt() {
  const d = EVT.res;
  const g = $("#evt-grandeur").value;
  const ecart = $("#evt-ecart-u").checked;
  const choisis = [...EVT.unites];
  const series = choisis.map((u, k) => {
    const y = d.series[g][u];
    const y0 = ecart ? y[0] ?? 0 : 0;
    return { xs: d.t, ys: y.map((v) => (v === null ? NaN : v - y0)), nom: u, couleur: TEINTES[k % TEINTES.length] };
  });
  tracerCourbes($("#evt-courbes-u"), "evtu", series, {
    titreX: "temps [s]", titreY: ecart ? `${g} - ${g}(0) [pu]` : `${g} [pu]`, verticales: verticalesEvt(d), zeroY: true,
  });
  const filtre = $("#evt-filtre-u").value.trim().toLowerCase();
  const l = d.unites.filter((u) => !filtre || u.nom.toLowerCase().includes(filtre));
  // les unites sorties d'abord, dans l'ordre ou elles sont sorties
  l.sort((a, b) => (a.sortie ?? Infinity) - (b.sortie ?? Infinity) || Math.abs(b.dQ ?? 0) - Math.abs(a.dQ ?? 0));
  tableEvt("#evt-table-unites",
    ["", "unite", "genre", "P0", "Q0", "P fin", "Q fin", "dP", "dQ", "sortie a [s]"],
    l.map((u) => {
      const k = choisis.indexOf(u.nom);
      const tr = el("tr", { class: [k >= 0 ? "choisi" : "", u.sortie !== null ? "sortie" : ""].join(" ") },
        el("td", {}, pastilleEvt(k)), el("td", { class: "mono" }, u.nom), el("td", {}, u.genre),
        celluleEvt(u.P0), celluleEvt(u.Q0), celluleEvt(u.Pfin), celluleEvt(u.Qfin),
        celluleEvt(u.dP), celluleEvt(u.dQ), celluleEvt(u.sortie, 3));
      tr.addEventListener("click", () => {
        if (EVT.unites.has(u.nom)) EVT.unites.delete(u.nom);
        else EVT.unites.add(u.nom);
        rendreUnitesEvt();
      });
      return tr;
    }));
}

function rendreFrequenceEvt() {
  const d = EVT.res;
  const f = d.frequence;
  const svg = $("#evt-courbe-f");
  if (!f.series.moy) {
    svg.textContent = "";
    $("#evt-f-titre").textContent = "pas de mesure";
    return;
  }
  const s = (cle) => f.series[cle].map((v) => v ?? NaN);
  tracerCourbes(svg, "evtf", [
    { xs: d.t, ys: s("min"), nom: "minimum des bus", couleur: "#1f4e79", largeur: 1.2 },
    { xs: d.t, ys: s("max"), nom: "maximum des bus", couleur: "#bf2617", largeur: 1.2 },
    { xs: d.t, ys: s("moy"), nom: "moyenne", couleur: "#1c2330", largeur: 2 },
  ], { titreX: "temps [s]", titreY: "f [Hz]", verticales: verticalesEvt(d) });
  const lesBus = f.bus.length > 8 ? `${f.bus.length} bus` : `bus ${f.bus.join(", ")}`;
  $("#evt-f-titre").textContent = `${lesBus} ; min ${fmt(f.min, 3)} Hz, max ${fmt(f.max, 3)} Hz, fin ${fmt(f.fin, 3)} Hz`;
}

function rendreRelaisEvt() {
  const d = EVT.res;
  const n = $("#evt-relais-liste");
  n.textContent = "";
  if (d.relais) {
    const ret = d.relais.retombee === null ? "aucune (verrouille)" : d.relais.retombee;
    n.append(el("p", { class: "doux" },
      `seuil ${d.relais.v_max} pu, temporisation ${d.relais.retard} s, retombee ${ret} ; ` +
      `${d.demarrages.length} demarrage(s), ${d.retombees.length} retombee(s)`));
  } else if (!d.declenchements.length) {
    n.append(el("p", { class: "doux" }, "aucun relais arme"));
    return;
  }
  n.append(el("div", { id: "evt-relais-table", class: "defile" }));
  tableEvt("#evt-relais-table", ["protection", "demarre a [s]", "declenche a [s]", "manoeuvres"],
    d.declenchements.map((x) => el("tr", {},
      el("td", { class: "mono" }, x.protection), celluleEvt(x.t_demarrage, 3), celluleEvt(x.t, 3),
      el("td", {}, x.manoeuvres.map((m) => `${m.genre} ${m.cible}`).join(", ")))));
}

async function cascadeEvt() {
  const etat = $("#evt-etat");
  const ouvertures = EVT.liste.filter((e) => ["ouverture", "blocage"].includes(e.genre));
  const initiales = [...new Set(ouvertures.map((e) => e.cible))];
  if (!initiales.length) {
    etat.textContent = "la cascade part d'ouvertures : ajouter l'ouverture (ou le blocage) d'au moins une unite";
    return;
  }
  const relais = relaisEvt() || { v_max: parseFloat($("#evt-vmax").value) || 1.1, familles: ["GFL"] };
  const v = E.version;
  etat.textContent = "cascade...";
  try {
    const d = await lancerTache("cascade", corpsReseau({
      initiales, v_max: relais.v_max, familles: relais.familles, ordre: $("#evt-ordre").value,
    }));
    if (v !== E.version) return;
    EVT.cascade = d;
    rendreCascadeEvt();
    const ignores = EVT.liste.length - ouvertures.length;
    etat.textContent = `cascade ${d.issue} : ${d.declenchees.length} unite(s) declenchee(s)` +
      (ignores ? ` ; ${ignores} evenement(s) autres que des ouvertures ignore(s)` : "");
  } catch (e) {
    etat.textContent = e.message;
  }
}

function rendreCascadeEvt() {
  const d = EVT.cascade;
  $("#evt-rangee-cascade").hidden = false;
  const xs = [0, ...d.etapes.map((e) => e.rang + 1)];
  const decl = new Set(d.declenchees.map((n) => d.bus_unite[n]));
  const init = new Set(d.initiales.map((n) => d.bus_unite[n]));
  const bus = [...new Set(Object.values(d.bus_unite))].sort((a, b) => a - b);
  const series = bus.map((b) => ({
    xs, ys: [d.V0[b] ?? NaN, ...d.etapes.map((e) => e.V[b] ?? NaN)],
    nom: `bus ${b}` + (decl.has(b) ? " (declenche)" : init.has(b) ? " (ouvert au depart)" : ""),
    couleur: decl.has(b) ? "#bf2617" : init.has(b) ? "#8a93a5" : "#1f4e79",
    largeur: decl.has(b) ? 1.8 : 1,
  }));
  series.push({ xs: [0, xs[xs.length - 1]], ys: [d.seuil, d.seuil], couleur: "#bf2617", pointille: true, largeur: 1.2, nom: `seuil ${d.seuil} pu` });
  tracerCourbes($("#evt-courbes-cascade"), "evtc", series, { titreX: "regime permanent (0 : depart)", titreY: "|V| [pu]" });
  const ordre = d.ordre === "vague" ? "en vagues" : "une a une";
  $("#evt-cascade-titre").textContent =
    `${ordre} ; ${d.issue} ; ${d.declenchees.length} declenchee(s)` + (d.declenchees.length ? ` : ${d.declenchees.join(", ")}` : "");
  const vmax0 = Math.max(...Object.values(d.V0));
  tableEvt("#evt-cascade-table", ["regime", "ouvertes a ce regime", "V max", "au bus", "au-dessus du seuil"], [
    el("tr", {}, el("td", {}, "0 (depart)"), el("td", {}, "--"), celluleEvt(vmax0), el("td", {}, "--"), el("td", {}, "--")),
    ...d.etapes.map((e) => el("tr", { class: e.converge ? "" : "instable" },
      el("td", {}, String(e.rang + 1)), el("td", { class: "mono" }, e.ouvertes.join(", ")),
      celluleEvt(e.v_max), el("td", { class: "mono" }, e.bus_v_max === null ? "--" : String(e.bus_v_max)),
      el("td", { class: "mono" }, e.converge ? (e.au_dessus.join(", ") || "aucune : la cascade s'arrete") : "load flow non converge"))),
  ]);
}

function activerOnglet(nom) {
  // les tables sont une VUE de l'onglet Reseau (creation par tables)
  const nav = nom === "tables" ? "constructeur" : nom;
  $$("nav button").forEach((b) => b.classList.toggle("actif", b.dataset.onglet === nav));
  $$(".bascule-creation button").forEach((b) => b.classList.toggle("actif", b.dataset.mode === nom));
  $$(".onglet").forEach((s) => s.classList.toggle("actif", s.id === "onglet-" + nom));
  if (nom === "constructeur") dessiner();
  if (nom === "tables") rendreTables();
  if (nom === "loadflow") {
    rendrePoint();
    if (E.res.lf) rendreLoadflow();
    if (E.res.scc) rendreScc();
  }
  if (nom === "impact" && E.res.impact) rendreImpact();
  if (nom === "robustesse" && E.res.robustesse) rendreRobustesse();
  if (nom === "simplification") {
    if (E.simpInvVersion !== E.version) chargerSimplification().then(() => E.res.simplification && rendreSimplification());
    else if (E.res.simplification) rendreSimplification();
  }
  if (nom === "modal" && E.res.modal) rendreModal();
  if (nom === "temporel" && !E.res.signaux) chargerSignaux().catch((e) => ($("#temporel-etat").textContent = e.message));
  if (nom === "evenements") chargerEvenements().catch((e) => ($("#evt-etat").textContent = e.message));
  if (nom === "balayage") {
    pret().catch(() => {}).then(rendreCibleBalayage);
    if (E.res.balayage && E.lieux) dessinerLieux();
  }
  if (nom === "blocs") rendreBlocs();
  if (nom === "equations") rendreEquations();
  if (nom === "doc" && !E.docDemandee) ouvrirDoc();
  if (nom === "stabilite") {
    remplirStatismes();
    chargerElements().catch((e) => ($("#cont-etat").textContent = e.message));
  }
}

async function rafraichirCatalogue(choisi) {
  const cat = await api("/api/reseaux");
  const sel = $("#catalogue");
  sel.textContent = "";
  for (const f of cat)
    sel.append(el("option", { value: f.nom }, `${f.nom} \u2014 ${f.bus} bus, ${f.sm} SM, ${f.gfm} GFM, ${f.gfl} GFL`));
  if (choisi) sel.value = choisi;
  return cat;
}

async function demarrer() {
  E.defauts = await api("/api/defauts");

  const pal = $("#palette-outils");
  for (const groupe of PALETTE) {
    pal.append(el("div", { class: "sep" }, groupe.titre));
    for (const o of groupe.outils) {
      const b = el("button", { "data-outil": o.outil, title: o.aide || o.titre }, iconePalette(o.icone, o.source));
      b.append(o.titre);
      pal.append(b);
    }
  }
  $$(".palette button[data-outil]").forEach((b) => b.addEventListener("click", () => outil(b.dataset.outil)));

  // la liste du catalogue arrive pendant qu'on branche les boutons : un clic
  // sur Nouveau ou Ouvrir des le premier affichage ne doit pas se perdre
  const catalogue = rafraichirCatalogue();
  catalogue.catch(() => {}); // l'erreur remonte au await de la fin
  const sel = $("#catalogue");

  $("#btn-charger").addEventListener("click", () => chargerCatalogue(sel.value));
  $("#btn-nouveau").addEventListener("click", () => {
    sortirZones();
    poser(reseauVide(), {}, "nouveau");
  });
  $("#btn-ouvrir").addEventListener("click", () => $("#fichier").click());
  $("#fichier").addEventListener("change", (e) => {
    if (e.target.files[0]) ouvrir(e.target.files[0]);
    e.target.value = ""; // rouvrir le meme fichier doit redeclencher
  });
  $("#btn-enregistrer").addEventListener("click", enregistrer);
  $("#btn-rapport").addEventListener("click", rapportEtude);
  $("#btn-catalogue").addEventListener("click", enregistrerAuCatalogue);
  $("#nom-reseau").addEventListener("change", (e) => {
    E.tables.nom = e.target.value;
    sauverBrouillon();
  });
  $("#btn-supprimer").addEventListener("click", supprimerSelection);
  $("#btn-copier").addEventListener("click", copier);
  $("#btn-import").addEventListener("click", ouvrirImport);
  $("#imp-source").addEventListener("change", changerSourceImport);
  $("#imp-fichier").addEventListener("change", changerSourceImport);
  $("#imp-mode").addEventListener("change", majModeImport);
  $("#btn-imp-choisir").addEventListener("click", () => outil("choix_import"));
  $("#btn-importer").addEventListener("click", importerReseau);
  $("#btn-imp-fermer").addEventListener("click", () => { $("#panneau-import").hidden = true; });
  $("#btn-coller").addEventListener("click", coller);
  $("#btn-annuler").addEventListener("click", () => revenir(-1));
  $("#btn-retablir").addEventListener("click", () => revenir(+1));
  document.addEventListener("keydown", (ev) => {
    const actif = document.activeElement || {};
    if (["INPUT", "TEXTAREA", "SELECT"].includes(actif.tagName) || actif.isContentEditable) return;
    if (!(ev.ctrlKey || ev.metaKey)) return;
    const k = ev.key.toLowerCase();
    if (k === "z" && !ev.shiftKey) {
      ev.preventDefault();
      revenir(-1);
    } else if (k === "y" || (k === "z" && ev.shiftKey)) {
      ev.preventDefault();
      revenir(+1);
    }
  });
  document.addEventListener("keydown", (ev) => {
    const ecrit = ["INPUT", "TEXTAREA", "SELECT"].includes((document.activeElement || {}).tagName);
    if (ecrit || !$("#onglet-constructeur").classList.contains("actif")) return;
    if ((ev.ctrlKey || ev.metaKey) && ev.key === "c") copier();
    else if ((ev.ctrlKey || ev.metaKey) && ev.key === "v") coller();
    else if (ev.key === "Escape") {
      viderGroupe();
      outil("select");
    }
  });
  $("#btn-placer").addEventListener("click", async () => {
    const d = await api("/api/topologie", { tables: E.tables, positions: null });
    E.positions = d.positions;
    dessiner();
  });
  $$("nav button").forEach((b) => b.addEventListener("click", () => activerOnglet(b.dataset.onglet)));
  $$(".bascule-creation button").forEach((b) => b.addEventListener("click", () => activerOnglet(b.dataset.mode)));
  $$('.bascule-creation button[data-mode="constructeur"]').forEach((b) => b.classList.add("actif"));
  $("#btn-loadflow").addEventListener("click", calculerLoadflow);
  $("#btn-scc").addEventListener("click", calculerScc);
  $("#btn-modal").addEventListener("click", calculerModal);
  $("#osc-seuls").addEventListener("change", rendreModal);
  $("#forme-suffixe").addEventListener("change", rendreForme);
  $("#mode-animer").addEventListener("change", rendreModeSchema);
  $("#filtre-domaine").addEventListener("change", rendreModal);
  $("#btn-impact-dc").addEventListener("click", calculerImpactContinu);
  $("#btn-stat").addEventListener("click", balayerStatismes);
  $("#btn-cont").addEventListener("click", suivreParametre);
  $("#btn-imp").addEventListener("click", boucleMineure);
  initFrise();
  $("#btn-rejeu").addEventListener("click", basculerRejeu);
  $("#rejeu-t").addEventListener("input", (e) => {
    arreterRejeu();
    rendreRejeu(+e.target.value);
  });
  $("#duree").addEventListener("change", rendreFrise);
  $("#imp-compresse").addEventListener("change", () => {
    ecrireBornes("nyq", null);
    rendreNyquist();
  });
  $("#cont-element").addEventListener("change", remplirParametresContinuation);
  $("#cont-param").addEventListener("change", bornesParametre);
  $("#btn-temporel").addEventListener("click", simuler);
  $("#evt-genre").addEventListener("change", () => remplirCiblesEvt(true));
  $("#evt-filtre").addEventListener("input", () => remplirCiblesEvt());
  $("#btn-evt-ajouter").addEventListener("click", ajouterEvt);
  $("#btn-evt-vider").addEventListener("click", () => {
    EVT.liste = [];
    rendreListeEvt();
  });
  $("#btn-evt-simuler").addEventListener("click", simulerEvt);
  $("#btn-evt-cascade").addEventListener("click", cascadeEvt);
  $("#evt-relais").addEventListener("change", (e) => ($("#evt-reglages-relais").hidden = !e.target.checked));
  for (const s of ["#evt-tri", "#evt-hors-seuls"]) $(s).addEventListener("change", () => EVT.res && rendreTensionsEvt());
  $("#evt-filtre-bus").addEventListener("input", () => EVT.res && rendreTensionsEvt());
  for (const s of ["#evt-grandeur", "#evt-ecart-u"]) $(s).addEventListener("change", () => EVT.res && rendreUnitesEvt());
  $("#evt-filtre-u").addEventListener("input", () => EVT.res && rendreUnitesEvt());
  $("#btn-balayage").addEventListener("click", balayer);
  $("#btn-robustesse").addEventListener("click", robustesse);
  $("#btn-simp-comparer").addEventListener("click", comparerSimplification);
  $("#btn-simp-effacer").addEventListener("click", () => {
    E.simp = {};
    simpChange();
  });
  $("#simp-appliquer").addEventListener("change", (ev) => {
    E.simpAppliquer = ev.target.checked;
    sauverBrouillon();
    E.version += 1;
    E.res = { lf: null, modal: null, signaux: null, forme: null, balayage: null, validation: null, impact: null };
    E.simpInvVersion = E.version;
  });
  $("#bal-curseur").addEventListener("input", (e) => {
    animerLieux(false);
    placerCurseurLieux(+e.target.value);
  });
  $("#bal-anim").addEventListener("click", () => animerLieux(!E._animLieux));
  $("#bal-element").addEventListener("change", () => remplirParametres(true));
  $("#bal-param").addEventListener("change", choisirParametre);
  $("#btn-ajouter-evt").addEventListener("click", ajouterEvenement);
  $("#btn-export").addEventListener("click", exporter);
  $$("[data-exporter]").forEach((b) => b.addEventListener("click", () => exporterVers(b.dataset.exporter)));
  $("#btn-code").addEventListener("click", produireCode);
  $("#btn-copier-code").addEventListener("click", () => navigator.clipboard && navigator.clipboard.writeText($("#code-python").value));
  $("#btn-impact").addEventListener("click", calculerImpact);
  $("#btn-mul-gfl").addEventListener("click", ajouterGflAleatoires);
  $("#btn-mul-comparer").addEventListener("click", comparerModelesReduits);
  $("#btn-penetration").addEventListener("click", balayerPenetration);
  $("#btn-doc-accueil").addEventListener("click", () => ouvrirDoc("index.html"));
  $("#btn-doc-equations").addEventListener("click", () => ouvrirDoc("equations_modeles.html"));
  $("#btn-doc-api").addEventListener("click", () => ouvrirDoc("api/index.html"));
  $("#btn-doc-construire").addEventListener("click", construireDoc);
  $("#btn-un-slack").addEventListener("click", unSlackToutPQ);
  $("#voir-transits").addEventListener("change", rendreLoadflow);
  for (const id of ["#axe-x", "#axe-y"]) {
    const sel = $(id);
    for (const [k, a] of Object.entries(AXES)) sel.append(el("option", { value: k }, a.nom));
    sel.addEventListener("change", () => choisirAxes($("#axe-x").value, $("#axe-y").value));
  }
  $("#axe-x").value = "re";
  $("#axe-y").value = "f";
  $$(".prereglage").forEach((b) => b.addEventListener("click", () => choisirAxes(...b.dataset.axes.split(","))));
  $("#btn-csv-teig").addEventListener("click", () =>
    E.res.modal && telecharger(csvTeig(modesVisibles(E.res.modal)), "text/csv", "T_eig.csv"));
  initApercuTables();
  document.addEventListener("keydown", (e) => {
    if (e.key === "Delete" && E.sel && document.activeElement === document.body) supprimerSelection();
  });

  initCanevas();
  surveillerFigures();
  outil("select");
  // le travail en cours d'abord ; le premier exemple du catalogue sinon --
  // sauf si l'on a deja pose un reseau pendant le demarrage (Nouveau, Ouvrir)
  const cat = await catalogue;
  if (!E.tables && !restaurerBrouillon()) await chargerCatalogue(cat.length ? cat[0].nom : "smib", true);
}

demarrer().catch((err) => etat("demarrage impossible : " + err.message, true));
