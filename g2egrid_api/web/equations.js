/* g2egrid -- les EQUATIONS dynamiques, modifiables (FR-052).

   Trois vues d'une meme chose, la liste d'equations ``dx/dt = f(x, u)`` :

   * l'onglet Equations : la liste GLOBALE du reseau, un equipement apres
     l'autre (filtres par equipement et par genre, recherche) ;
   * la section Equations de la fiche et des proprietes d'un equipement : la
     meme liste restreinte a l'equipement, avec « Ajouter un etat »,
     « Ajouter une variable » et « Revenir au modele natif » ;
   * le panneau d'un modele a equations dans l'onglet Blocs.

   Ce fichier ne calcule rien : le serveur (/api/equations/...) convertit,
   verifie et ecrit. Toute modification revient dans E.tables par
   modifie() : annuler / retablir la couvrent, et la validation du reseau
   repart apres chaque ecriture. Aucune expression ni aucun nom n'entre dans
   la page par innerHTML : tout passe par textContent. */

"use strict";

const ROUTES_EQ = {
  reseau: "/api/equations/reseau",
  element: "/api/equations/element",
  modifier: "/api/equations/modifier",
  ajouter: "/api/equations/ajouter",
  supprimer: "/api/equations/supprimer",
  parametre: "/api/equations/parametre",
  revenir: "/api/equations/revenir",
  verifier: "/api/equations/verifier",
  // passages entre l'onglet Blocs et les equations
  depuisSchema: "/api/equations/depuis_schema",
  versBlocs: "/api/equations/vers_blocs",
};

const Q = {
  lignes: null, // la liste globale (reponse de /api/equations/reseau)
  versionLignes: -1, // E.version de ce calcul
  elem: "", // filtre par equipement ("" : tous)
  genres: new Set(["etat"]), // filtre par genre : les etats seuls (dx/dt = f) par defaut
  recherche: "",
  latex: false, // affichage LaTeX plutot que texte
  limite: 300, // lignes affichees (une liste de grand reseau est longue)
  ouvert: { fiche: true, inspecteur: false }, // section Equations repliee ou non
  params: { fiche: false, inspecteur: false }, // parametres du modele : replies ou non, par section
  paramsModele: false, // de meme pour le modele a equations de l'onglet Blocs
  sections: new Set(), // les sections vivantes (fiche, proprietes), a redessiner
  chargement: 0,
};

/** Les tables dont un equipement a des equations modifiables (port alternatif). */
const TABLES_EQUATIONS = /^T_(SMs|GFM\ds|GFL\ds|IBs|Loads|Blocs)$/;

const NOMS_RESERVES = ["wb", "kS", "vd", "vq", "V", "w", "theta_g"];

const GENRES_EQ = {
  etat: ["etat", "variable d'etat : la liste donne dx/dt = f(x, u)"],
  algebrique: ["algebrique", "variable algebrique : z = g(x, u) (ou 0 = g si elle est implicite)"],
  sortie: ["sortie", "sortie vers le reseau : courant injecte id, iq (pu base reseau), ou P, Q ; angle theta et vitesse w d'un formeur"],
};

// ----------------------------------------------------------------- appels
/** Un appel POST aux routes /api/equations : comme ``appelBlocs``, mais les
 *  messages disent quoi corriger (route absente, expression refusee). */
async function appelEquations(chemin, corps) {
  let r;
  try {
    r = await fetch(chemin, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(corps),
    });
  } catch {
    throw new Error("le serveur ne repond pas : verifier qu'il tourne toujours");
  }
  const txt = await r.text();
  let data = null;
  try {
    data = JSON.parse(txt);
  } catch {
    // pas du JSON : le statut parlera
  }
  if (r.ok) return data;
  const detail = data && data.detail;
  const msg = typeof detail === "string"
    ? detail
    : Array.isArray(detail)
      ? detail.map((d) => `${(d.loc || []).slice(1).join(".") || "requete"} : ${d.msg}`).join(" ; ")
      : txt.slice(0, 200);
  if (r.status === 404 && (!detail || detail === "Not Found"))
    throw new Error(`les equations modifiables n'existent pas dans la version du serveur lancee (route ${chemin} absente) : relancer l'atelier avec la version a jour`);
  if (r.status === 422) throw new Error(msg || "le serveur refuse cette equation");
  if (r.status === 403) throw new Error("acces refuse : rouvrir l'atelier par son lien avec le jeton");
  throw new Error(msg || `erreur ${r.status}`);
}

/** Les tables rendues par le serveur remplacent celles du reseau (comme une
 *  reparation) : les modeles et les lignes de T_Blocs en viennent. */
function appliquerTablesEquations(nouvelles) {
  const nom = E.tables.nom;
  const t = { ...nouvelles };
  delete t.SIMPLIFICATIONS; // celles de l'onglet restent les seules
  // le serveur n'ecrit pas un MODELES_BLOCS vide : sans cela, le modele d'un
  // equipement rendu a son etat natif resterait dans le reseau
  if (!t.MODELES_BLOCS) delete E.tables.MODELES_BLOCS;
  Object.assign(E.tables, t);
  E.tables.nom = nom;
}

// ------------------------------------------------- validation apres ecriture
/** Le reseau tient-il encore ? (controles de structure, puis load flow, dont
 *  l'initialisation des equipements). Le message reste affiche tant qu'une
 *  autre ecriture ne l'a pas remplace. */
function afficherValidationEquations(etatV, texte) {
  for (const n of $$(".eq-validation")) {
    n.dataset.etat = etatV;
    n.className = "eq-validation pb " + (etatV === "ok" ? "ok" : etatV === "attente" ? "avis" : "erreur");
    n.textContent = texte;
    n.hidden = !texte;
  }
}

async function validerApresEquation() {
  const v = E.version;
  afficherValidationEquations("attente", "validation du reseau...");
  try {
    const d = await api("/api/valider", corpsReseau());
    if (v !== E.version) return;
    const erreurs = d.problemes.filter((p) => p.gravite === "erreur");
    if (erreurs.length) {
      const m = "Validation : " + erreurs.map((p) => p.message).join(" ; ");
      afficherValidationEquations("erreur", m);
      etat(m, true);
      return;
    }
    const lf = await api("/api/loadflow", corpsReseau());
    if (v !== E.version) return;
    if (lf.converge === false) {
      const m = "Load flow non converge avec ces equations : verifier les consignes et les valeurs de depart des variables (annuler : Ctrl+Z)";
      afficherValidationEquations("erreur", m);
      etat(m, true);
      return;
    }
    afficherValidationEquations("ok", "Reseau valide, load flow converge.");
  } catch (err) {
    if (v !== E.version) return;
    const m = "Validation : " + err.message;
    afficherValidationEquations("erreur", m);
    etat(m, true);
  }
}

// ------------------------------------- confirmation : fonctionnement libre
const TEXTE_GARDES =
  "Le modele genere decrit le fonctionnement libre : la limitation de courant et la saturation " +
  "de tension du modele natif (gardes) ne sont pas reprises. L'element ne se limitera plus en " +
  "simulation temporelle. L'operation reste annulable (Ctrl+Z).";

/** Demande confirmation avant de remplacer un equipement natif par un modele
 *  (blocs ou equations) : rend une promesse (true : on continue). Boite
 *  construite par el(), sans innerHTML ; Echap et le fond annulent. */
function confirmerFonctionnementLibre(quoi) {
  return new Promise((resolve) => {
    const fin = (ok) => {
      fond.remove();
      document.removeEventListener("keydown", touche, true);
      resolve(ok);
    };
    const touche = (ev) => {
      if (ev.key === "Escape") {
        ev.stopPropagation();
        ev.preventDefault();
        fin(false);
      }
    };
    const annuler = el("button", { type: "button", "data-action": "annuler-gardes" }, "Annuler");
    const ok = el("button", { type: "button", class: "principal", "data-action": "confirmer-gardes" }, "Remplacer quand meme");
    annuler.addEventListener("click", () => fin(false));
    ok.addEventListener("click", () => fin(true));
    const fond = el("div", { class: "confirm-fond", "data-confirmation": "gardes" },
      el("div", { class: "confirm-boite", role: "alertdialog", "aria-modal": "true", "aria-label": "Fonctionnement libre" },
        el("h3", {}, "Les limitations ne sont pas reprises"),
        el("p", {}, quoi),
        el("p", { class: "confirm-texte" }, TEXTE_GARDES),
        el("div", { class: "confirm-actions" }, annuler, " ", ok)));
    fond.addEventListener("mousedown", (ev) => {
      if (ev.target === fond) fin(false);
    });
    document.addEventListener("keydown", touche, true);
    document.body.append(fond);
    annuler.focus();
  });
}

// ---------------------------------------------------------------- ecriture
/** Les tables rendues par une ecriture remplacent celles du reseau. Rend la
 *  reference de l'equipement dans ces tables : le serveur la donne
 *  (``table``, ``index``) -- un natif converti passe dans T_Blocs, un
 *  equipement rendu a son modele natif retrouve sa table. */
function appliquerEcriture(d, ref) {
  if (!d || !d.tables) throw new Error("reponse du serveur sans tables");
  appliquerTablesEquations(d.tables);
  const nouvelle = typeof d.table === "string" && Number.isInteger(d.index)
    ? { table: d.table, index: d.index }
    : ref;
  suivreReference(ref, nouvelle);
  modifie();
  validerApresEquation();
  return nouvelle;
}

/** Une ecriture (modifier, ajouter, supprimer, parametre) sur l'equipement
 *  ``ref`` : le serveur rend les tables et le modele. Rend la nouvelle
 *  reference de l'equipement (un natif qui devient modele a equations change
 *  de table). Lance l'erreur du serveur telle quelle, les tables intactes. */
async function ecrireEquation(op, ref, extra) {
  if (!E.tables) throw new Error("charger ou creer un reseau d'abord");
  if (!lignes(ref.table)[ref.index]) throw new Error(`${ref.table} n'a pas de ligne ${ref.index}`);
  if (ref.table !== "T_Blocs" && !(await confirmerFonctionnementLibre(
    `${ref.table}[${ref.index}] sera remplace par son modele a equations pour appliquer cette modification.`)))
    throw new Error("modification annulee : l'equipement reste natif");
  const d = await appelEquations(ROUTES_EQ[op], { tables: E.tables, table: ref.table, index: ref.index, ...extra });
  return appliquerEcriture(d, ref);
}

/** L'equipement a change de table : la selection et la fiche le suivent. */
function suivreReference(ancienne, nouvelle) {
  if (ancienne.table === nouvelle.table && ancienne.index === nouvelle.index) return;
  if (E.sel && E.sel.table === ancienne.table && E.sel.index === ancienne.index) E.sel = { ...nouvelle };
  if (E.ficheRef && E.ficheRef.table === ancienne.table && E.ficheRef.index === ancienne.index)
    E.ficheRef = { ...nouvelle };
}

/** Rend l'equipement de T_Blocs[index] a sa table d'origine. C'est le
 *  SERVEUR qui le fait : le modele garde son origine (la ligne native, la
 *  puissance d'une charge reportee a son bus), y compris dans le fichier du
 *  reseau -- le retour marche donc apres un rechargement, et pour une charge. */
async function revenirAuNatif(ref) {
  const ligneB = lignes("T_Blocs")[ref.index];
  if (!ligneB) throw new Error("cet equipement n'est plus dans le reseau");
  const nom = ligneB[2];
  const d = await appelEquations(ROUTES_EQ.revenir, { tables: E.tables, table: "T_Blocs", index: ref.index });
  const nouvelle = appliquerEcriture(d, ref);
  if (typeof B !== "undefined" && B.modele === nom && !modelesBlocs()[nom]) B.modele = null;
  return nouvelle;
}

// -------------------------------------------------------------------- LaTeX
function nomLatex(n) {
  const m = /^([A-Za-z]+)_(.+)$/.exec(n);
  return m ? `${m[1]}_{${m[2].replace(/_/g, ",")}}` : n;
}

/** Un apercu LaTeX d'une expression ecrite avec les noms du modele (le
 *  serveur peut fournir le sien dans ``latex``). */
function expressionLatex(e) {
  return String(e)
    .replace(/\*\*/g, "^")
    .replace(/\*/g, " \\cdot ")
    .replace(/\b([A-Za-z]+)_([A-Za-z0-9_]+)/g, (_, a, b) => `${a}_{${b.replace(/_/g, ",")}}`);
}

function gaucheTexte(l) {
  return l.genre === "etat" ? `d${l.variable}/dt` : l.implicite ? "0" : l.variable;
}

function gaucheLatex(l) {
  const v = nomLatex(l.variable);
  return l.genre === "etat" ? `\\frac{d ${v}}{dt}` : l.implicite ? "0" : v;
}

function droiteLatex(l) {
  return l.latex ? String(l.latex).replace(/^.*?=\s*/, "") : expressionLatex(l.expression);
}

function listeLatex(rows) {
  return rows.map((l) => `${gaucheLatex(l)} = ${droiteLatex(l)}`).join(" \\\\\n");
}

function copierLatex(rows, noeud) {
  const t = listeLatex(rows);
  try {
    if (navigator.clipboard) navigator.clipboard.writeText(t).catch(() => {});
  } catch {
    // presse-papiers refuse : le texte reste visible ci-dessous
  }
  if (noeud) {
    noeud.textContent = "";
    noeud.append(el("pre", { class: "mono eq-latex-copie" }, t));
  }
}

// ------------------------------------------------------------ la liste (table)
const CAPTION_EXPR = "Expression f(x, u) : les noms des etats, des variables algebriques, des consignes et des parametres du modele, et les reserves wb, kS (base) et les entrees du reseau vd, vq, V, w, theta_g. Grandeurs en pu.";

/** Le tableau d'equations. ``opts`` : ``modifier(ligne, expression)`` (promesse,
 *  lance l'erreur a montrer), ``supprimer(ligne)`` (facultatif), ``avecElement``
 *  (colonnes equipement et bus), ``noms`` (liste d'autocompletion). */
function tableEquations(rows, opts = {}) {
  const t = el("table", { class: "eq-table" });
  const tete = el("tr", {});
  if (opts.avecElement) {
    tete.append(el("th", { title: "Equipement du reseau" }, "Equipement"), el("th", { title: "Bus de l'equipement" }, "Bus"));
  }
  tete.append(
    el("th", { title: "A gauche la variable : pour un etat, sa derivee dx/dt", class: "eq-g" }, "Variable"),
    el("th", { class: "eq-egal" }, ""),
    el("th", { class: "eq-d", title: CAPTION_EXPR }, "Expression (clic pour modifier)"),
  );
  if (opts.supprimer) tete.append(el("th", {}, ""));
  t.append(el("thead", {}, tete));
  const corps = el("tbody", {});
  for (const l of rows) corps.append(ligneEquation(l, opts));
  t.append(corps);
  return t;
}

function ligneEquation(l, opts) {
  const g = GENRES_EQ[l.genre] || [l.genre, ""];
  const tr = el("tr", {
    class: "eq-ligne" + (l.modifiable ? "" : " eq-grise") + (l.modifiee ? " eq-modifiee" : ""),
    "data-element": l.element || "",
    "data-variable": l.variable,
    "data-genre": l.genre,
    "data-modifiable": l.modifiable ? "1" : "0",
    "data-modifiee": l.modifiee ? "1" : "0",
    title: l.modifiable ? "" : (l.pourquoi || "equation non modifiable"),
  });
  if (opts.avecElement) {
    const nom = el("td", { class: "eq-nom" }, el("span", {}, l.element || ""));
    if (l.modifiee) nom.append(" ", el("span", { class: "eq-badge", title: "equipement dont le modele a ete remplace par un modele a equations" }, "modifie"));
    tr.append(nom, el("td", { class: "mono" }, l.bus === null || l.bus === undefined ? "" : String(l.bus)));
  }
  tr.append(
    el("td", { class: "eq-var mono eq-g", title: g[1] }, Q.latex ? gaucheLatex(l) : gaucheTexte(l)),
    el("td", { class: "eq-egal", title: l.implicite ? "variable implicite : 0 = expression" : "" }, "="),
  );
  const td = el("td", { class: "eq-expr mono eq-d" });
  affichageExpression(td, l);
  if (l.modifiable && opts.modifier) {
    td.tabIndex = 0;
    td.title = "Cliquer pour modifier (Entree valide, Echap annule)";
    td.addEventListener("click", () => ouvrirEdition(td, l, opts));
    td.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter" && ev.target === td) ev.preventDefault(), ouvrirEdition(td, l, opts);
    });
  }
  tr.append(td);
  if (opts.supprimer) {
    const sup = el("td", {});
    if (l.genre !== "sortie" && l.modifiable)
      sup.append(el("button", {
        type: "button", class: "petit eq-supprimer", "data-supprimer": l.variable,
        title: `Supprimer ${l.variable} (refuse si une autre equation s'en sert)`,
        onclick: () => opts.supprimer(l),
      }, "✕"));
    tr.append(sup);
  }
  return tr;
}

function affichageExpression(td, l) {
  td.textContent = Q.latex ? droiteLatex(l) : l.expression;
  if (l.implicite) td.dataset.implicite = "1";
}

/** L'edition sur place : un champ a la place de l'expression. Entree valide
 *  (le serveur verifie avant d'ecrire), Echap annule ; une erreur reste
 *  affichee sous le champ, qui garde ce qu'on a tape. */
function ouvrirEdition(td, l, opts) {
  if (td.querySelector("input")) return;
  const champ = el("input", {
    class: "eq-champ mono", value: l.expression, spellcheck: "false", autocomplete: "off",
    "aria-label": `equation de ${l.variable}`,
  });
  if (opts.noms) champ.setAttribute("list", opts.noms);
  const msg = el("div", { class: "eq-erreur pb erreur", hidden: "" });
  td.textContent = "";
  td.append(champ, msg);
  let occupe = false;
  const fermer = () => {
    if (!td.isConnected) return;
    td.textContent = "";
    affichageExpression(td, l);
    td.focus({ preventScroll: true });
  };
  champ.addEventListener("keydown", async (ev) => {
    if (ev.key === "Escape") {
      ev.stopPropagation();
      ev.preventDefault();
      fermer();
    } else if (ev.key === "Enter") {
      ev.preventDefault();
      if (occupe) return;
      const v = champ.value.trim();
      if (!v) {
        msg.textContent = "l'expression est vide : ecrire une expression, ou Echap pour annuler";
        msg.hidden = false;
        return;
      }
      if (v === l.expression) return fermer();
      occupe = true;
      champ.disabled = true;
      msg.hidden = true;
      try {
        await opts.modifier(l, v);
      } catch (err) {
        occupe = false;
        if (!td.isConnected) return;
        champ.disabled = false;
        msg.textContent = `${l.variable} : ${err.message}`;
        msg.hidden = false;
        champ.focus();
      }
    }
  });
  champ.addEventListener("blur", () => {
    setTimeout(() => {
      if (!occupe && !champ.disabled && td.isConnected && msg.hidden) fermer();
    }, 120);
  });
  champ.focus();
  champ.select();
}

// ----------------------------------------------------------- l'onglet global
function elementsDeLaListe() {
  const vus = new Map();
  for (const l of Q.lignes || []) if (!vus.has(l.element)) vus.set(l.element, l);
  return [...vus.values()];
}

function lignesFiltrees() {
  const rech = Q.recherche.trim().toLowerCase();
  return (Q.lignes || []).filter((l) => {
    if (!Q.genres.has(l.genre)) return false;
    if (Q.elem && l.element !== Q.elem) return false;
    if (!rech) return true;
    return `${l.element} ${l.variable} ${l.expression}`.toLowerCase().includes(rech);
  });
}

async function chargerEquationsReseau(force = false) {
  if (!E.tables) return false;
  if (!force && Q.lignes && Q.versionLignes === E.version) return true;
  const n = ++Q.chargement;
  $("#eq-etat").textContent = "calcul des equations...";
  $("#eq-etat").className = "doux";
  for (let essai = 0; essai < 3; essai++) {
    const v = E.version;
    try {
      const d = await appelEquations(ROUTES_EQ.reseau, { tables: E.tables });
      if (n !== Q.chargement) return false; // une demande plus recente a pris le relais
      if (v !== E.version) continue; // le reseau a change pendant le calcul
      Q.lignes = d.lignes || [];
      Q.versionLignes = v;
      $("#eq-etat").textContent = "";
      return true;
    } catch (err) {
      if (n !== Q.chargement) return false;
      Q.lignes = null;
      Q.versionLignes = -1;
      $("#eq-etat").textContent = err.message;
      $("#eq-etat").className = "pb erreur";
      return false;
    }
  }
  return false;
}

async function rendreEquations(force = false) {
  const racine = $("#onglet-equations");
  if (!racine) return;
  if (!E.tables) {
    $("#eq-liste").textContent = "";
    $("#eq-compte").textContent = "";
    return;
  }
  if (force || !Q.lignes || Q.versionLignes !== E.version) {
    const ok = await chargerEquationsReseau(force);
    if (!ok && !Q.lignes) {
      $("#eq-liste").textContent = "";
      $("#eq-compte").textContent = "";
      return;
    }
  }
  dessinerEquations();
}

function dessinerEquations() {
  // le filtre par equipement
  const sel = $("#eq-element");
  const elems = elementsDeLaListe();
  if (Q.elem && !elems.some((l) => l.element === Q.elem)) Q.elem = "";
  sel.textContent = "";
  sel.append(el("option", { value: "" }, `tous (${elems.length})`));
  for (const l of elems)
    sel.append(el("option", { value: l.element, selected: l.element === Q.elem ? "" : null },
      `${l.modifiee ? "● " : ""}${l.element} (${l.classe || l.table || ""})`));
  for (const g of Object.keys(GENRES_EQ)) $(`#eq-g-${g}`).checked = Q.genres.has(g);
  $("#eq-latex").checked = Q.latex;
  const toutes = Q.lignes || [];
  const rows = lignesFiltrees();
  const n = $("#eq-compte");
  const nbModif = elems.filter((l) => l.modifiee).length;
  n.textContent = `${rows.length} equation(s) sur ${toutes.length}` + (nbModif ? ` · ${nbModif} equipement(s) modifie(s)` : "");
  const liste = $("#eq-liste");
  liste.textContent = "";
  if (!rows.length) {
    liste.append(el("p", { class: "doux" }, toutes.length
      ? "Aucune equation ne correspond aux filtres (genres cochees, equipement, recherche)."
      : "Aucune equation : le reseau n'a pas d'equipement assemble."));
    return;
  }
  const vues = rows.slice(0, Q.limite);
  liste.append(tableEquations(vues, {
    avecElement: true,
    modifier: (l, expr) => modifierDepuisListe(l, expr),
    noms: "eq-noms-reseau",
  }));
  if (rows.length > vues.length) {
    const plus = el("button", { type: "button", id: "eq-plus" }, `Afficher ${Math.min(300, rows.length - vues.length)} equations de plus (${rows.length - vues.length} masquees)`);
    plus.addEventListener("click", () => {
      Q.limite += 300;
      dessinerEquations();
    });
    liste.append(plus);
  }
  remplirNomsReseau(rows);
}

/** L'autocompletion du champ d'edition : variables de la liste, reserves. */
function remplirNomsReseau(rows) {
  let dl = $("#eq-noms-reseau");
  if (!dl) {
    dl = el("datalist", { id: "eq-noms-reseau" });
    document.body.append(dl);
  }
  dl.textContent = "";
  const noms = new Set(NOMS_RESERVES);
  for (const l of rows) noms.add(l.variable);
  for (const n of noms) dl.append(el("option", { value: n }));
}

async function modifierDepuisListe(l, expr) {
  const nouvelle = await ecrireEquation("modifier", { table: l.table, index: l.index }, { variable: l.variable, expression: expr });
  etat(`equation de ${l.variable} (${l.element}) modifiee : annulable par Ctrl+Z`);
  await rafraichirEquations(true);
  return nouvelle;
}

/** Tout ce qui montre des equations se met a jour. */
async function rafraichirEquations(force = false) {
  if ($("#onglet-equations").classList.contains("actif")) await rendreEquations(force);
  else if (force) Q.lignes = null;
  for (const s of [...Q.sections]) {
    if (!s.det.isConnected) Q.sections.delete(s);
    else await s.rafraichir();
  }
}

function installerEquations() {
  for (const g of Object.keys(GENRES_EQ))
    $(`#eq-g-${g}`).addEventListener("change", (ev) => {
      if (ev.target.checked) Q.genres.add(g);
      else Q.genres.delete(g);
      Q.limite = 300;
      dessinerEquations();
    });
  $("#eq-element").addEventListener("change", (ev) => {
    Q.elem = ev.target.value;
    Q.limite = 300;
    dessinerEquations();
  });
  let minuteur = null;
  $("#eq-recherche").addEventListener("input", (ev) => {
    clearTimeout(minuteur);
    minuteur = setTimeout(() => {
      Q.recherche = ev.target.value;
      Q.limite = 300;
      if (Q.lignes) dessinerEquations();
    }, 120);
  });
  $("#eq-latex").addEventListener("change", (ev) => {
    Q.latex = ev.target.checked;
    if (Q.lignes) dessinerEquations();
    for (const s of Q.sections) s.dessiner();
    if (typeof rendreModeleEquations === "function") rendreModeleEquations();
  });
  $("#eq-copier").addEventListener("click", () => copierLatex(lignesFiltrees(), $("#eq-latex-sortie")));
  $("#eq-actualiser").addEventListener("click", () => rendreEquations(true));
  $("#btn-blocs-convertir").addEventListener("click", convertirModeleBlocs);
}

// ----------------------------------------- section d'un equipement (fiche, proprietes)
/** La section Equations de l'equipement ``table[index]`` : la meme liste,
 *  restreinte a lui. ``cle`` : "fiche" ou "inspecteur" (l'etat replie/ouvert
 *  se retient par cle). */
function sectionEquations(table, index, cle) {
  const det = el("details", {
    class: "avance eq-section", "data-section-equations": cle, open: Q.ouvert[cle] ? "" : null,
  }, el("summary", { title: "Les equations dynamiques dx/dt = f(x, u) de cet equipement ; modifiables sur place" },
    "Equations modifiables (dx/dt = f(x, u))"));
  const corps = el("div", { class: "eq-corps" });
  det.append(corps);
  const ctl = { det, ref: { table, index }, modele: null, info: {}, charge: -1, erreur: null };

  ctl.rafraichir = async () => {
    if (!det.open || !det.isConnected) return;
    const v = E.version;
    corps.replaceChildren(el("p", { class: "doux" }, "lecture des equations..."));
    try {
      const d = await appelEquations(ROUTES_EQ.element, { tables: E.tables, table: ctl.ref.table, index: ctl.ref.index });
      if (!det.isConnected) return;
      ctl.modele = d.modele;
      ctl.info = { natif: d.natif, schemaBloc: !!d.schema_bloc };
      ctl.erreur = null;
      ctl.charge = v;
    } catch (err) {
      if (!det.isConnected) return;
      ctl.modele = null;
      ctl.erreur = err.message;
    }
    ctl.dessiner();
  };

  ctl.dessiner = () => {
    corps.textContent = "";
    if (ctl.erreur) {
      corps.append(el("div", { class: "pb avis eq-pas-modifiable" }, "Equations non modifiables ou indisponibles : " + ctl.erreur));
      return;
    }
    const m = ctl.modele;
    if (!m) return;
    const natif = ctl.ref.table !== "T_Blocs";
    // ``natif`` rendu par le serveur : le nom de l'equipement natif remplace,
    // null si ce modele n'en remplace aucun (ecrit a la main : pas de retour)
    const origine = ctl.info.natif;
    const rows = lignesDeModele(m, ctl.ref, true);
    const nbEt = rows.filter((l) => l.genre === "etat").length;
    const quoi = ctl.info.schemaBloc ? "Schema-bloc" : "Modele a equations";
    corps.append(el("p", { class: "doux eq-resume" },
      natif
        ? `Modele natif ${m.nom || ""} ecrit en equations : ${nbEt} etat(s). La premiere modification remplace l'equipement par ce modele a equations (meme bus, meme puissance) ; « Revenir au modele natif » le rend.`
        : `${quoi} ${m.nom || ""} : ${nbEt} etat(s), ${rows.length} equation(s). ` +
          (origine === null
            ? "Ce modele ne remplace aucun equipement natif."
            : `Il remplace l'equipement natif${origine ? " " + origine : ""} ; « Revenir au modele natif » le rend.`) +
          (ctl.info.schemaBloc ? " Modifier une equation en fait un modele a equations a part : le dessin reste dans l'onglet Blocs." : "")));
    corps.append(tableEquations(rows, {
      modifier: (l, expr) => ctl.modifier(l, expr),
      supprimer: (l) => ctl.supprimer(l),
      noms: "eq-noms-" + cle,
    }));
    const dl = el("datalist", { id: "eq-noms-" + cle });
    for (const n of nomsDuModele(m)) dl.append(el("option", { value: n }));
    corps.append(dl);
    const barre = el("div", { class: "eq-barre" });
    barre.append(
      el("button", { type: "button", "data-action": "ajouter-etat", title: "Ajoute un etat dx/dt = f(x, u), avec sa valeur de depart" }, "Ajouter un etat"),
      " ",
      el("button", { type: "button", "data-action": "ajouter-variable", title: "Ajoute une variable algebrique z = g(x, u), avec sa valeur de depart" }, "Ajouter une variable"),
    );
    for (const b of barre.querySelectorAll("button[data-action^=ajouter]"))
      b.addEventListener("click", () => ctl.formulaireAjout(b.dataset.action === "ajouter-etat" ? "etat" : "algebrique"));
    if (!natif && origine !== null) {
      const rev = el("button", {
        type: "button", "data-action": "revenir-natif",
        title: "L'equipement retrouve sa table d'origine : le modele a equations est retire du reseau",
      }, "Revenir au modele natif");
      rev.addEventListener("click", () => ctl.revenir());
      barre.append(" ", rev);
    }
    const cp = el("button", { type: "button", "data-action": "copier-latex" }, "Copier en LaTeX");
    const sortieLatex = el("div", {});
    cp.addEventListener("click", () => copierLatex(rows, sortieLatex));
    const lat = el("input", { type: "checkbox", "data-action": "latex" });
    lat.checked = Q.latex;
    lat.addEventListener("change", () => {
      Q.latex = lat.checked;
      $("#eq-latex").checked = Q.latex;
      for (const s of Q.sections) s.dessiner();
      if (Q.lignes && $("#onglet-equations").classList.contains("actif")) dessinerEquations();
    });
    barre.append(" ", cp, " ", el("label", { class: "doux" }, lat, " LaTeX"));
    corps.append(barre);
    corps.append(el("div", { class: "eq-ajout", hidden: "" }));
    const par = blocParametres(m.parametres, natif
      ? { ouvert: Q.params[cle], note: "Equipement natif : ses parametres se reglent dans ses proprietes (« Reglages »), il reste alors natif. Ici, en lecture seule." }
      : { ouvert: Q.params[cle], modifier: (nom, valeur) => ctl.parametre(nom, valeur), note: "Une valeur validee (Entree) est ecrite dans le modele de cet equipement ; annulable par Ctrl+Z." });
    if (par) {
      par.addEventListener("toggle", () => {
        Q.params[cle] = par.open;
      });
      corps.append(par);
    }
    corps.append(el("div", { class: "eq-validation pb", hidden: "" }), sortieLatex);
  };

  ctl.modifier = async (l, expr) => {
    ctl.ref = await ecrireEquation("modifier", ctl.ref, { variable: l.variable, expression: expr });
    etat(`equation de ${l.variable} modifiee : annulable par Ctrl+Z`);
    await rafraichirEquations(true);
  };

  ctl.supprimer = async (l) => {
    try {
      ctl.ref = await ecrireEquation("supprimer", ctl.ref, { variable: l.variable });
      etat(`${l.variable} supprime : annulable par Ctrl+Z`);
      await rafraichirEquations(true);
    } catch (err) {
      montrerErreurSection(corps, `Supprimer ${l.variable} : ${err.message}`);
    }
  };

  ctl.parametre = async (nom, valeur) => {
    ctl.ref = await ecrireEquation("parametre", ctl.ref, { nom, valeur });
    etat(`parametre ${nom} = ${valeur} : annulable par Ctrl+Z`);
    await rafraichirEquations(true);
  };

  ctl.revenir = async () => {
    try {
      ctl.ref = await revenirAuNatif(ctl.ref);
      etat("equipement rendu a son modele natif : annulable par Ctrl+Z");
      await rafraichirEquations(true);
    } catch (err) {
      montrerErreurSection(corps, "Retour au modele natif impossible : " + err.message);
    }
  };

  ctl.formulaireAjout = (genre) => {
    const zone = corps.querySelector(".eq-ajout");
    if (!zone) return;
    zone.hidden = false;
    zone.textContent = "";
    const etatG = genre === "etat";
    const nom = el("input", { class: "mono", size: 12, placeholder: etatG ? "nom de l'etat" : "nom de la variable", "data-champ": "nom", "aria-label": "nom" });
    const expr = el("input", { class: "mono eq-champ", placeholder: etatG ? "derivee, ex. wb*(w - 1)" : "expression, ex. 2*V", "data-champ": "expression", "aria-label": "expression", list: "eq-noms-" + cle });
    const dep = el("input", { type: "number", step: "any", value: "0", style: "width:7em", "data-champ": "depart", title: "Valeur de depart du load flow (point d'equilibre initial)", "aria-label": "valeur de depart" });
    const ok = el("button", { type: "button", class: "primaire", "data-action": "valider-ajout" }, "Ajouter");
    const annuler = el("button", { type: "button" }, "Annuler");
    const msg = el("div", { class: "pb erreur", hidden: "" });
    const valider = async () => {
      const v = nom.value.trim();
      if (!/^[A-Za-z][A-Za-z0-9_]*$/.test(v)) {
        msg.textContent = "nom invalide : lettres, chiffres et _, commencant par une lettre";
        msg.hidden = false;
        return;
      }
      if (!expr.value.trim()) {
        msg.textContent = etatG ? "ecrire la derivee de l'etat" : "ecrire l'expression de la variable";
        msg.hidden = false;
        return;
      }
      const extra = { genre, variable: v, expression: expr.value.trim() };
      const d0 = parseFloat(dep.value);
      if (!Number.isNaN(d0)) extra.depart = d0;
      ok.disabled = true;
      try {
        ctl.ref = await ecrireEquation("ajouter", ctl.ref, extra);
        etat(`${etatG ? "etat" : "variable"} ${v} ajoute : annulable par Ctrl+Z`);
        await rafraichirEquations(true);
      } catch (err) {
        ok.disabled = false;
        msg.textContent = err.message;
        msg.hidden = false;
      }
    };
    ok.addEventListener("click", valider);
    annuler.addEventListener("click", () => ctl.dessiner());
    zone.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") ev.preventDefault(), valider();
      else if (ev.key === "Escape") ev.stopPropagation(), ctl.dessiner();
    });
    zone.append(
      el("b", {}, etatG ? "Nouvel etat : " : "Nouvelle variable : "),
      nom, etatG ? " d/dt = " : " = ", expr, " depart ", dep, " ", ok, " ", annuler, msg);
    nom.focus();
  };

  det.addEventListener("toggle", () => {
    Q.ouvert[cle] = det.open;
    if (det.open && ctl.charge !== E.version) ctl.rafraichir();
  });
  Q.sections.add(ctl);
  if (det.open) ctl.rafraichir();
  return det;
}

function montrerErreurSection(corps, texte) {
  for (const n of corps.querySelectorAll(".eq-erreur-section")) n.remove();
  corps.append(el("div", { class: "pb erreur eq-erreur-section" }, texte));
}

/** Les parametres d'un modele, nom et valeur, dans une section repliee.
 *  ``opts.modifier(nom, valeur)`` (promesse) rend les valeurs modifiables :
 *  une valeur refusee revient, avec la raison du serveur. Sans lui, la liste
 *  est en lecture seule. Rend null s'il n'y a aucun parametre. */
function blocParametres(parametres, opts = {}) {
  const noms = Object.keys(parametres || {});
  if (!noms.length) return null;
  const det = el("details", { class: "eq-parametres", "data-section-parametres": "", open: opts.ouvert ? "" : null },
    el("summary", { title: "Les valeurs numeriques des parametres que citent les equations" }, `Parametres du modele (${noms.length})`));
  if (opts.note) det.append(el("p", { class: "doux" }, opts.note));
  const grille = el("div", { class: "eq-params-grille" });
  // le message d'un refus n'existe que le temps du refus (pas de bloc cache)
  const effacer = () => {
    for (const n of det.querySelectorAll(".eq-param-erreur")) n.remove();
  };
  for (const n of noms) {
    const v = parametres[n];
    if (!opts.modifier) {
      grille.append(el("span", { class: "mono eq-param", "data-parametre": n }, `${n} = ${fmt(v, 6)}`));
      continue;
    }
    const champ = el("input", { type: "number", step: "any", value: String(v), class: "mono", "data-parametre": n, "aria-label": `parametre ${n}` });
    champ.addEventListener("change", async () => {
      const x = parseFloat(champ.value);
      effacer();
      if (!Number.isFinite(x)) {
        champ.value = String(v);
        return;
      }
      if (x === +v) return;
      champ.disabled = true;
      try {
        await opts.modifier(n, x);
      } catch (err) {
        champ.disabled = false;
        champ.value = String(v);
        det.append(el("div", { class: "pb erreur eq-param-erreur" }, `${n} : ${err.message}`));
      }
    });
    grille.append(el("label", { class: "eq-param" }, el("span", { class: "mono" }, n), " = ", champ));
  }
  det.append(grille);
  return det;
}

/** Les equations d'un modele a equations (``{etats, algebriques, sorties}``)
 *  sous la forme d'une ligne par equation, comme la liste globale. */
function lignesDeModele(m, ref, modifiable) {
  const base = { table: ref.table, index: ref.index, modifiable, modifiee: ref.table === "T_Blocs", element: m.nom || "" };
  const rows = [];
  for (const e of m.etats || []) rows.push({ ...base, genre: "etat", variable: e.nom, expression: String(e.derivee) });
  for (const a of m.algebriques || []) {
    const implicite = !!(a.implicite || (a.residu !== undefined && a.expr === undefined));
    rows.push({ ...base, genre: "algebrique", variable: a.nom, expression: String(implicite ? a.residu : a.expr), implicite });
  }
  for (const [k, v] of Object.entries(m.sorties || {})) rows.push({ ...base, genre: "sortie", variable: k, expression: String(v) });
  return rows;
}

function nomsDuModele(m) {
  const noms = new Set(NOMS_RESERVES);
  for (const e of m.etats || []) noms.add(e.nom);
  for (const a of m.algebriques || []) noms.add(a.nom);
  for (const c of m.consignes || []) noms.add(typeof c === "string" ? c : c.nom);
  for (const p of Object.keys(m.parametres || {})) noms.add(p);
  return noms;
}

// ---------------------------------- un modele a equations dans l'onglet Blocs
/** Le panneau qui remplace le dessin quand le modele ouvert est un modele a
 *  equations : ses equations (modifiables, verifiees par le serveur avant
 *  d'etre gardees) et le passage vers les blocs. */
function rendreModeleEquations() {
  const s = typeof schemaCourant === "function" ? schemaCourant() : null;
  const p = $("#blocs-eq-modele");
  const btn = $("#btn-blocs-convertir");
  const eq = !!s && s.genre === "equations";
  p.hidden = !eq;
  btn.textContent = eq ? "Convertir en blocs" : "Voir en equations";
  btn.title = eq
    ? "Cree un schema-bloc equivalent (nouveau modele) et l'ouvre"
    : "Cree le modele a equations equivalent (nouveau modele) et l'ouvre : variable a gauche, expression a droite";
  if (!eq) return;
  p.textContent = "";
  const rows = lignesDeModele(s, { table: "MODELES_BLOCS", index: 0 }, true).map((l) => ({ ...l, modifiee: false }));
  p.append(el("h3", {}, `Modele a equations ${s.nom}`));
  p.append(el("p", { class: "doux" },
    "Les equations dx/dt = f(x, u) de ce modele. Une modification est verifiee (analyseur sur, noms connus, compilation) avant d'etre gardee. " +
    "Les equipements qui utilisent ce modele la suivent."));
  p.append(tableEquations(rows, { modifier: (l, expr) => modifierModeleEquations(s, l, expr), noms: "eq-noms-modele" }));
  const par = blocParametres(s.parametres, {
    ouvert: Q.paramsModele,
    modifier: (nom, valeur) => modifierParametreModele(s, nom, valeur),
    note: "Une valeur validee (Entree) est verifiee par une compilation avant d'etre gardee.",
  });
  if (par) {
    par.addEventListener("toggle", () => {
      Q.paramsModele = par.open;
    });
    p.append(par);
  }
  const dl = el("datalist", { id: "eq-noms-modele" });
  for (const n of nomsDuModele(s)) dl.append(el("option", { value: n }));
  p.append(dl);
  const cp = el("button", { type: "button", "data-action": "copier-latex" }, "Copier en LaTeX");
  const sortieLatex = el("div", {});
  cp.addEventListener("click", () => copierLatex(rows, sortieLatex));
  p.append(el("div", { class: "eq-barre" }, cp), el("div", { class: "eq-validation pb", hidden: "" }), sortieLatex);
}

/** Corriger une equation dans le modele lui-meme : le serveur compile une
 *  copie a blanc (/verifier) avant que le modele ne change. */
async function modifierModeleEquations(s, l, expr) {
  const copie = structuredClone(s);
  let trouve = false;
  if (l.genre === "etat") {
    const e = (copie.etats || []).find((x) => x.nom === l.variable);
    if (e) e.derivee = expr, trouve = true;
  } else if (l.genre === "algebrique") {
    const a = (copie.algebriques || []).find((x) => x.nom === l.variable);
    if (a) {
      if (a.implicite || (a.residu !== undefined && a.expr === undefined)) a.residu = expr;
      else a.expr = expr;
      trouve = true;
    }
  } else if (copie.sorties && l.variable in copie.sorties) {
    copie.sorties[l.variable] = expr;
    trouve = true;
  }
  if (!trouve) throw new Error(`${l.variable} n'est plus dans ce modele`);
  const d = await appelEquations(ROUTES_EQ.verifier, { modele: copie });
  if (d && d.ok === false) throw new Error(d.message || "modele refuse par la compilation");
  modelesBlocs()[s.nom] = copie;
  modifie();
  rendreModeleEquations();
  etat(`equation de ${l.variable} modifiee dans le modele ${s.nom} : annulable par Ctrl+Z`);
}

/** La valeur d'un parametre dans le modele lui-meme (onglet Blocs) : comme
 *  une equation, une copie est compilee a blanc avant que le modele change. */
async function modifierParametreModele(s, nom, valeur) {
  const copie = structuredClone(s);
  if (!copie.parametres || !(nom in copie.parametres)) throw new Error(`${nom} n'est plus un parametre de ce modele`);
  copie.parametres[nom] = valeur;
  const d = await appelEquations(ROUTES_EQ.verifier, { modele: copie });
  if (d && d.ok === false) throw new Error(d.message || "modele refuse par la compilation");
  modelesBlocs()[s.nom] = copie;
  modifie();
  rendreModeleEquations();
  etat(`parametre ${nom} = ${valeur} dans le modele ${s.nom} : annulable par Ctrl+Z`);
}

/** « Voir en equations » (schema-bloc) ou « Convertir en blocs » (modele a
 *  equations) : un NOUVEAU modele, le modele d'origine n'est pas touche. */
async function convertirModeleBlocs() {
  const s = typeof schemaCourant === "function" ? schemaCourant() : null;
  if (!s) return etat("aucun modele ouvert", true);
  const versBlocs = s.genre === "equations";
  try {
    const d = versBlocs
      ? await appelEquations(ROUTES_EQ.versBlocs, { modele: s })
      : await appelEquations(ROUTES_EQ.depuisSchema, { schema_bloc: s });
    const nouveau = versBlocs ? d.schema : d.modele;
    if (!nouveau) throw new Error("reponse du serveur sans modele");
    nouveau.nom = nomLibre(`${s.nom}_${versBlocs ? "blocs" : "eq"}`);
    modelesBlocs()[nouveau.nom] = nouveau;
    ouvrirModeleBlocs(nouveau.nom);
    modifie();
    await rendreBlocs();
    etat(`modele ${nouveau.nom} cree depuis ${s.nom}`);
  } catch (err) {
    const m = (versBlocs ? "Conversion en blocs : " : "Passage en equations : ") + err.message;
    etat(m, true);
    $("#blocs-etat").textContent = m;
  }
}

/** « Compiler » d'un modele a equations : la compilation a blanc du serveur. */
async function verifierModeleEquations() {
  const s = schemaCourant();
  const sortie = $("#blocs-eq-modele");
  for (const n of sortie.querySelectorAll(".eq-verif")) n.remove();
  try {
    const d = await appelEquations(ROUTES_EQ.verifier, { modele: s });
    const ok = d.ok !== false;
    sortie.append(el("div", { class: "pb eq-verif " + (ok ? "ok" : "erreur") },
      ok ? `Modele compile : ${(d.equations || []).length} equation(s).` : d.message || "modele refuse"));
  } catch (err) {
    sortie.append(el("div", { class: "pb erreur eq-verif" }, err.message));
  }
}

installerEquations();
