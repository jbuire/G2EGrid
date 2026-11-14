// =====================================================================
// Onglet Blocs : un petit Simulink pour decrire un nouvel element du reseau
// =====================================================================
// Un MODELE est un schema-bloc (g2egrid.blocs) : des blocs (consigne, gain,
// somme, 1er ordre, PI...) relies par des liens, deux sorties vers le reseau
// (P et Q, ou id et iq). Les modeles vivent dans E.tables.MODELES_BLOCS et
// voyagent donc avec le reseau (fichier, annuler/retablir) ; un element de
// T_Blocs designe son modele par son nom. Le serveur compile le schema en
// equations (/api/blocs/compiler) : c'est ce que le moteur assemblera.
//
// La position d'un bloc est rangee dans le bloc lui-meme (x, y) : le
// compilateur ignore ces champs.

const B = {
  modele: null, // nom du modele ouvert
  sel: null, // {genre: "bloc" | "lien", id | index}
  types: null, // /api/blocs/types : types de blocs et exemples
  pret: null, // la promesse de leur chargement
  drag: null, // deplacement d'un bloc ou trace d'un lien en cours
  compile: null, // dernier resultat de compilation
  chemin: [], // identifiants des sous-systemes ouverts : le niveau affiche
  multi: new Set(), // les blocs selectionnes (du niveau affiche)
  clic: null, // le dernier clic sur un bloc (detection du double-clic)
};
const BLOC_L = 92, BLOC_H = 50;

/** Les exemples viennent du serveur (B.types.exemples, g2egrid.blocs) : un
 *  seul endroit les decrit. Celui-ci ne sert que si l'outil de la palette du
 *  reseau pose un element avant que l'onglet ait jamais repondu. */
const SOURCE_DE_REPLI = {
  nom: "source_filtree",
  parametres: { tau: 0.05 },
  blocs: [
    { id: "pr", type: "consigne", nom: "P_ref", x: 40, y: 60 },
    { id: "qr", type: "consigne", nom: "Q_ref", x: 40, y: 200 },
    { id: "fp", type: "premier_ordre", K: 1, tau: "tau", x: 230, y: 60 },
    { id: "fq", type: "premier_ordre", K: 1, tau: "tau", x: 230, y: 200 },
    { id: "op", type: "sortie", signal: "P", x: 430, y: 60 },
    { id: "oq", type: "sortie", signal: "Q", x: 430, y: 200 },
  ],
  liens: [
    { de: "pr", vers: "fp", entree: 0 }, { de: "qr", vers: "fq", entree: 0 },
    { de: "fp", vers: "op", entree: 0 }, { de: "fq", vers: "oq", entree: 0 },
  ],
};

const COULEUR_BLOC = {
  entree: "#00743f", sortie: "#1f4e79", consigne: "#bf7c17", constante: "#667085",
};

function modelesBlocs() {
  if (!E.tables) return {};
  return E.tables.MODELES_BLOCS || (E.tables.MODELES_BLOCS = {});
}

/** Le modele que l'outil de la palette du reseau pose : celui ouvert dans
 *  l'onglet, sinon le premier, sinon un exemple cree pour l'occasion. */
function modeleBlocsCourant() {
  const m = modelesBlocs();
  if (B.modele && m[B.modele]) return B.modele;
  const noms = Object.keys(m);
  if (noms.length) return (B.modele = noms[0]);
  const s = structuredClone(B.types ? B.types.exemples.source_filtree.schema : SOURCE_DE_REPLI);
  m[s.nom] = s;
  etat(`modele ${s.nom} cree (onglet Blocs pour le modifier)`);
  return (B.modele = s.nom);
}

/** Le modele ``nom`` ecrit-il son angle (sortie theta) ? Il peut alors
 *  porter le repere du reseau : c'est un formeur. */
function modeleBlocsFormeur(nom) {
  const s = modelesBlocs()[nom];
  if (s && s.genre === "equations") return !!(s.sorties && "theta" in s.sorties);
  return !!s && (s.blocs || []).some((b) => b.type === "sortie" && b.signal === "theta");
}

function schemaCourant() {
  const m = modelesBlocs();
  return B.modele ? m[B.modele] : null;
}

/** Les types propres a l'interface : ils completent ceux du serveur (qui les
 *  decrit aussi, avec leur categorie) tant qu'il ne les envoie pas, et fixent
 *  le titre de repli d'un type inconnu. */
const TYPES_LOCAUX = {
  sous_systeme: { titre: "Sous-systeme", categorie: "structure", champs: ["nom"], entrees: null, sorties: [], noms_entrees: [] },
  port_entree: { titre: "Port d'entree", categorie: "structure", champs: ["index", "nom"], entrees: 0, sorties: ["y"], noms_entrees: [] },
  port_sortie: { titre: "Port de sortie", categorie: "structure", champs: ["index", "nom"], entrees: 1, sorties: [], noms_entrees: [] },
};

/** Les types qui n'ont de sens qu'a l'interieur d'un sous-systeme. */
const TYPES_INTERNES = ["port_entree", "port_sortie"];

function specType(type) {
  return (B.types && B.types.types[type]) || TYPES_LOCAUX[type] || null;
}

/** Le niveau affiche : le schema du modele, ou celui du sous-systeme dans
 *  lequel on est entre (meme forme : blocs et liens). Un chemin devenu
 *  caduc (modele change, annulation) est ramene au dernier niveau valide. */
function niveau() {
  let n = schemaCourant();
  if (!n) return null;
  const valide = [];
  for (const id of B.chemin) {
    const b = (n.blocs || []).find((x) => x.id === id && x.type === "sous_systeme" && x.schema);
    if (!b) break;
    n = b.schema;
    valide.push(id);
  }
  B.chemin = valide;
  n.blocs = n.blocs || [];
  n.liens = n.liens || [];
  return n;
}

/** Le bloc sous-systeme ouvert (le dernier du chemin), ou null. */
function sousSystemeOuvert() {
  let n = schemaCourant(), b = null;
  for (const id of B.chemin) {
    b = ((n && n.blocs) || []).find((x) => x.id === id);
    n = b && b.schema;
  }
  return b;
}

/** Les ports d'un sous-systeme, dans l'ordre de leur ``index`` : c'est
 *  l'interface qu'il montre (la regle est celle du serveur ; ceci ne fait
 *  que la lire pour dessiner). */
function portsDe(b, type) {
  const blocs = (b.schema && b.schema.blocs) || [];
  return blocs
    .filter((x) => x.type === type)
    .sort((p, q) => (+p.index || 0) - (+q.index || 0))
    .map((x) => String(x.nom || x.id));
}

function nbEntrees(b) {
  if (b.type === "sous_systeme") return portsDe(b, "port_entree").length;
  const t = specType(b.type);
  if (!t) return 0;
  if (t.entrees !== null) return t.entrees;
  if (b.type === "somme") return String(b.signes || "++").length;
  return Math.max(1, +b.n || 1);
}

/** Les noms des sorties d'un bloc : une seule ("y") pour la plupart, "d"
 *  et "q" pour les blocs du repere dq, aucune pour une sortie vers le reseau,
 *  ceux des ports de sortie pour un sous-systeme. */
function nomsSorties(b) {
  if (b.type === "sous_systeme") return portsDe(b, "port_sortie");
  const t = specType(b.type);
  return t ? t.sorties : ["y"];
}

/** Les noms des entrees (affiches a cote des ports). */
function nomsEntrees(b) {
  if (b.type === "sous_systeme") return portsDe(b, "port_entree");
  const t = specType(b.type);
  return (t && t.noms_entrees) || [];
}
const nbSorties = (b) => nomsSorties(b).length;
const aSortie = (b) => nbSorties(b) > 0;
// une somme est un cercle : sa boite est carree
// les blocs du repere dq ont des noms de ports des deux cotes : plus larges
const largeurBloc = (b) => (b.type === "somme" ? BLOC_H : b.type === "sous_systeme" ? 130 : ["rl_dq", "c_dq", "park", "park_inverse"].includes(b.type) ? 120 : BLOC_L);
// un bloc a beaucoup de ports s'allonge
const hauteurBloc = (b) => (b.type === "somme" ? BLOC_H : Math.max(BLOC_H, 15 * (Math.max(nbEntrees(b), nbSorties(b)) + 1)));

function portEntree(b, k) {
  const n = nbEntrees(b), h = hauteurBloc(b);
  const y = b.y + (h * (k + 1)) / (n + 1);
  if (b.type !== "somme") return [b.x, y];
  // sur le cercle, cote gauche
  const r = BLOC_H / 2 - 2, cy = b.y + BLOC_H / 2;
  return [b.x + BLOC_H / 2 - Math.sqrt(Math.max(0, r * r - (y - cy) ** 2)), y];
}

function portSortie(b, k = 0) {
  const n = Math.max(1, nbSorties(b)), h = hauteurBloc(b);
  return [b.x + largeurBloc(b) - (b.type === "somme" ? 2 : 0), b.y + (h * (k + 1)) / (n + 1)];
}

/** Le texte d'un bloc, a la Simulink : la fonction de transfert. */
function libelleBloc(b) {
  const v = (x, d = "") => (x === undefined || x === "" ? d : String(x));
  switch (b.type) {
    case "entree": return v(b.signal, "?");
    case "sortie": return v(b.signal, "?");
    case "consigne": return v(b.nom, b.id);
    case "constante": return v(b.valeur, "0");
    case "gain": return v(b.K, "K");
    case "somme": return "Σ";
    case "produit": return "×";
    case "division": return "÷";
    case "fonction": return "f: " + v(b.expr, "u");
    case "integrateur": return `${v(b.K, "1")}/s`;
    case "premier_ordre": return `${v(b.K, "1")}/(1+${v(b.tau, "τ")}·s)`;
    case "pi": return `${v(b.Kp, "Kp")}+${v(b.Ki, "Ki")}/s`;
    case "derivee": return `${v(b.K, "1")}·s/(1+${v(b.tau, "τ")}·s)`;
    case "park": return "Park";
    case "park_inverse": return "Park⁻¹";
    case "rl_dq": return `1/(${v(b.R, "R")}+${v(b.L, "L")}·s)`;
    case "c_dq": return `1/(${v(b.C, "C")}·s)`;
    case "limiteur_courant": return `|i| ≤ ${v(b.Imax, "Imax")}`;
    case "saturation_tension": return `|m| ≤ ${v(b.Mmax, "Mmax")}`;
    case "pi_aw": return `${v(b.Kp, "Kp")}+${v(b.Ki, "Ki")}/s ⟲`;
    case "pid": return `PID ${v(b.Kp, "Kp")}+${v(b.Ki, "Ki")}/s+${v(b.Kd, "Kd")}·s`;
    case "avance_retard": return `${v(b.K, "K")}(1+${v(b.T1, "T1")}·s)/(1+${v(b.T2, "T2")}·s)`;
    case "second_ordre": return `${v(b.K, "K")}·${v(b.wn, "ωn")}²/(s²+2${v(b.zeta, "ζ")}${v(b.wn, "ωn")}s+${v(b.wn, "ωn")}²)`;
    case "fonction_transfert": return `${polynomeS(b.num)} / ${polynomeS(b.den)}`;
    case "saturation": return `[${v(b.min, "min")} ; ${v(b.max, "max")}]`;
    case "implicite": return "0 = " + v(b.expr, "f(u, y)");
    case "sous_systeme": return v(b.nom, b.id);
    case "port_entree": return `In ${(+b.index || 0) + 1} ${v(b.nom, "")}`.trim();
    case "port_sortie": return `Out ${(+b.index || 0) + 1} ${v(b.nom, "")}`.trim();
    default: return b.type;
  }
}

/** Un polynome en s, coefficients par puissances decroissantes (nombres ou
 *  expressions des parametres) : [1, "a", 0] -> "s²+a·s". */
function polynomeS(coefs) {
  const c = Array.isArray(coefs) ? coefs : [];
  if (!c.length) return "?";
  const sup = ["", "", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"];
  const termes = [];
  c.forEach((a, k) => {
    const p = c.length - 1 - k;
    const n = typeof a === "number" || /^[+-]?\d/.test(String(a)) ? a : null;
    if (n !== null && +n === 0) return;
    const txt = String(a);
    const simple = /^[A-Za-z0-9_.]+$/.test(txt) || n !== null;
    const coef = p > 0 && n !== null && +n === 1 ? "" : simple ? txt : `(${txt})`;
    const s = p === 0 ? "" : p === 1 ? "s" : "s" + (sup[p] ?? `^${p}`);
    termes.push(p === 0 || coef === "" ? coef + s : `${coef}·${s}`);
  });
  const t = termes.join("+").replace(/\+-/g, "-") || "0";
  return t.includes("+") || t.includes("-") ? `(${t})` : t;
}

// ------------------------------------------------------------------ rendu
/** Les types de blocs et les exemples, demandes une fois au serveur. */
function chargerTypesBlocs() {
  if (!B.pret) {
    B.pret = api("/api/blocs/types").then((d) => {
      B.types = d;
      rendrePaletteBlocs();
    });
    B.pret.catch(() => (B.pret = null));
  }
  return B.pret;
}

async function rendreBlocs() {
  if (!E.tables) return;
  try {
    await chargerTypesBlocs();
  } catch (err) {
    $("#blocs-etat").textContent = err.message;
    return;
  }
  const m = modelesBlocs();
  if (!B.modele || !m[B.modele]) B.modele = Object.keys(m)[0] || null;
  const sel = $("#blocs-modele");
  sel.textContent = "";
  for (const nom of Object.keys(m)) sel.append(el("option", { value: nom, selected: nom === B.modele ? "" : null }, nom));
  $("#blocs-vide").hidden = !!B.modele;
  const parEquations = !!B.modele && m[B.modele].genre === "equations";
  $("#blocs-espace").hidden = !B.modele || parEquations;
  if (typeof rendreModeleEquations === "function") rendreModeleEquations();
  if (!B.modele || parEquations) return;
  rendrePaletteBlocs();
  rendreAriane();
  dessinerSchemaBlocs();
  rendreProprietesBloc();
  rendreParametresModele();
  majBoutonsGroupe();
}

/** Les categories de la palette, dans l'ordre d'affichage : leur identifiant
 *  est celui que le serveur donne a chaque type (``categorie``). */
const CATEGORIES = [
  ["sources", "Sources et consignes"],
  ["operations", "Operations"],
  ["dynamique", "Dynamiques"],
  ["repere_dq", "Repere dq et filtre"],
  ["limiteurs", "Limitation"],
  ["structure", "Structure"],
  ["sorties", "Sorties vers le reseau"],
];

/** Repli tant que le serveur ne donne pas la categorie d'un type : le
 *  classement historique de la palette. Un type absent d'ici va dans "autres". */
const CATEGORIE_REPLI = {
  entree: "sources", consigne: "sources", constante: "sources", sortie: "sorties",
  gain: "operations", somme: "operations", produit: "operations", division: "operations", fonction: "operations",
  integrateur: "dynamique", premier_ordre: "dynamique", pi: "dynamique", derivee: "dynamique",
  park: "repere_dq", park_inverse: "repere_dq", rl_dq: "repere_dq", c_dq: "repere_dq",
  limiteur_courant: "limiteurs", saturation_tension: "limiteurs", pi_aw: "limiteurs",
};

const AIDE_TYPE = {
  entree: "lit une grandeur du reseau au bus de l'element (tension, frequence)",
  sortie: "ecrit une grandeur vers le reseau (P et Q, ou id et iq)",
  consigne: "entree de service, calee au load flow",
  pid: "Kp + Ki/s + Kd.s/(1+Tf.s)",
  avance_retard: "K.(1+T1.s)/(1+T2.s)",
  second_ordre: "K.wn²/(s²+2.zeta.wn.s+wn²)",
  fonction_transfert: "num(s)/den(s), coefficients par puissances decroissantes",
  saturation: "borne la sortie entre min et max (libre au point de fonctionnement)",
  implicite: "y tel que expr(u1..un, y) = 0 : contrainte algebrique",
  sous_systeme: "bloc qui en contient d'autres ; double-clic pour y entrer",
  port_entree: "entree du sous-systeme ouvert (son index fixe l'ordre des entrees)",
  port_sortie: "sortie du sous-systeme ouvert",
};

function rendrePaletteBlocs() {
  const p = $("#blocs-palette");
  p.textContent = "";
  if (!B.types) return;
  const dedans = B.chemin.length > 0;
  // tous les types connus, serveur d'abord, puis ceux que l'interface ajoute
  const tous = { ...TYPES_LOCAUX, ...B.types.types };
  const parCat = {};
  for (const [t, spec] of Object.entries(tous)) {
    if (TYPES_INTERNES.includes(t) && !dedans) continue;
    const cat = spec.categorie || CATEGORIE_REPLI[t] || "autres";
    (parCat[cat] = parCat[cat] || []).push(t);
  }
  const ordre = [...CATEGORIES, ...Object.keys(parCat).filter((c) => !CATEGORIES.some(([k]) => k === c)).map((c) => [c, c === "autres" ? "Autres" : c])];
  for (const [cat, titre] of ordre) {
    if (!parCat[cat]) continue;
    p.append(el("div", { class: "blocs-groupe", "data-categorie": cat }, titre));
    for (const t of parCat[cat]) {
      const spec = tous[t];
      const bt = el("button", { type: "button", "data-bloc": t, title: AIDE_TYPE[t] ? `${spec.titre} : ${AIDE_TYPE[t]}` : spec.titre }, spec.titre);
      bt.addEventListener("click", () => ajouterBloc(t));
      p.append(bt);
    }
  }
}

/** Le fil d'Ariane : modele > sous-systeme > ... ; un clic remonte. */
function rendreAriane() {
  const n = $("#blocs-ariane");
  if (!n) return;
  n.textContent = "";
  const s = schemaCourant();
  if (!s) return;
  const pas = [{ nom: s.nom || B.modele, k: 0 }];
  let niv = s;
  B.chemin.forEach((id, i) => {
    const b = niv.blocs.find((x) => x.id === id);
    pas.push({ nom: (b && (b.nom || b.id)) || id, k: i + 1 });
    niv = b && b.schema ? b.schema : { blocs: [] };
  });
  pas.forEach((a, i) => {
    if (i) n.append(el("span", { class: "doux", "aria-hidden": "true" }, " > "));
    const dernier = i === pas.length - 1;
    const bt = el("button", {
      type: "button", class: "ariane-pas" + (dernier ? " actif" : ""), "data-niveau": a.k,
      title: dernier ? "niveau affiche" : "remonter a ce niveau",
    }, String(a.nom));
    if (!dernier) bt.addEventListener("click", () => allerAuNiveau(a.k));
    n.append(bt);
  });
}

/** Changer de niveau : la selection n'a plus de sens, la palette change
 *  (les ports n'existent que dans un sous-systeme). */
function allerAuNiveau(k) {
  B.chemin = B.chemin.slice(0, k);
  B.sel = null;
  B.multi = new Set();
  B.clic = null;
  rendrePaletteBlocs();
  rendreAriane();
  dessinerSchemaBlocs();
  rendreProprietesBloc();
  majBoutonsGroupe();
}

function entrerSousSysteme(id) {
  const niv = niveau();
  const b = niv && niv.blocs.find((x) => x.id === id);
  if (!b || b.type !== "sous_systeme") return;
  if (!b.schema) b.schema = { blocs: [], liens: [] };
  B.chemin = [...B.chemin, id];
  B.sel = null;
  B.multi = new Set();
  B.clic = null;
  rendrePaletteBlocs();
  rendreAriane();
  dessinerSchemaBlocs();
  rendreProprietesBloc();
  majBoutonsGroupe();
}

function dessinerSchemaBlocs() {
  const s = niveau();
  const svg = $("#blocs-svg");
  svg.textContent = "";
  if (!s) return;
  const defs = el("g");
  svg.append(defs);
  const parId = Object.fromEntries(s.blocs.map((b) => [b.id, b]));
  // les liens d'abord, sous les blocs
  s.liens.forEach((l, i) => {
    const a = parId[l.de], c = parId[l.vers];
    if (!a || !c) return;
    const [x1, y1] = portSortie(a, +l.sortie || 0), [x2, y2] = portEntree(c, +l.entree || 0);
    const xm = x2 > x1 + 20 ? (x1 + x2) / 2 : x1 + 15;
    const d = x2 > x1 + 20
      ? `M ${x1} ${y1} H ${xm} V ${y2} H ${x2}`
      : `M ${x1} ${y1} H ${xm} V ${Math.max(a.y + hauteurBloc(a), c.y + hauteurBloc(c)) + 18} H ${x2 - 15} V ${y2} H ${x2}`;
    const choisi = B.sel && B.sel.genre === "lien" && B.sel.index === i;
    const g = el("g", { class: "bloc-lien", "data-lien": i });
    g.append(el("path", { d, fill: "none", stroke: "transparent", "stroke-width": 10 }));
    g.append(el("path", { d, fill: "none", stroke: choisi ? ROUGE_SEL : "#344054", "stroke-width": choisi ? 2.4 : 1.6 }));
    g.append(el("polygon", { points: `${x2},${y2} ${x2 - 7},${y2 - 4} ${x2 - 7},${y2 + 4}`, fill: choisi ? ROUGE_SEL : "#344054" }));
    g.addEventListener("mousedown", (ev) => {
      ev.stopPropagation();
      B.sel = { genre: "lien", index: i };
      B.multi = new Set();
      dessinerSchemaBlocs();
      rendreProprietesBloc();
      majBoutonsGroupe();
    });
    svg.append(g);
  });
  for (const b of s.blocs) svg.append(dessinerBloc(b));
  // la vue suit le schema : un GFL complet deborde largement 1000 x 560
  if (!B.drag) {
    const W = Math.max(1000, ...s.blocs.map((b) => b.x + largeurBloc(b) + 60));
    const H = Math.max(560, ...s.blocs.map((b) => b.y + hauteurBloc(b) + 50));
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  }
  if (B.drag && B.drag.genre === "lien" && B.drag.fin) {
    const [x1, y1] = portSortie(parId[B.drag.de], B.drag.sortie || 0);
    const [x2, y2] = B.drag.fin;
    svg.append(el("line", { x1, y1, x2, y2, stroke: ROUGE_SEL, "stroke-width": 1.5, "stroke-dasharray": "4 3" }));
  }
  if (B.drag && B.drag.genre === "rect" && B.drag.fin) {
    const [x0, y0] = B.drag.debut, [x1, y1] = B.drag.fin;
    svg.append(el("rect", {
      class: "blocs-rectangle", x: Math.min(x0, x1), y: Math.min(y0, y1), width: Math.abs(x1 - x0), height: Math.abs(y1 - y0),
      fill: "rgba(31,78,121,0.08)", stroke: "#1f4e79", "stroke-dasharray": "4 3", "pointer-events": "none",
    }));
  }
}

function dessinerBloc(b) {
  const choisi = B.multi.has(b.id);
  const sous = b.type === "sous_systeme";
  const encre = choisi ? ROUGE_SEL : sous ? "#5b3a8c" : COULEUR_BLOC[b.type] || "#344054";
  const g = el("g", { class: "bloc" + (sous ? " bloc-sous-systeme" : "") + (choisi ? " bloc-choisi" : ""), "data-bloc-id": b.id });
  const rond = ["entree", "sortie", "consigne", "constante", "port_entree", "port_sortie"].includes(b.type);
  const h = hauteurBloc(b);
  if (b.type === "gain") {
    g.append(el("polygon", {
      points: `${b.x},${b.y} ${b.x + BLOC_L},${b.y + BLOC_H / 2} ${b.x},${b.y + BLOC_H}`,
      fill: "#fff", stroke: encre, "stroke-width": choisi ? 2.4 : 1.6,
    }));
  } else if (b.type === "somme") {
    g.append(el("rect", { x: b.x, y: b.y, width: BLOC_H, height: BLOC_H, fill: "transparent" }));
    g.append(el("circle", { cx: b.x + BLOC_H / 2, cy: b.y + BLOC_H / 2, r: BLOC_H / 2 - 2, fill: "#fff", stroke: encre, "stroke-width": choisi ? 2.4 : 1.6 }));
  } else {
    if (sous) {
      // un sous-systeme est un bloc double : on le reconnait d'un coup d'oeil
      g.append(el("rect", { x: b.x + 4, y: b.y + 4, width: largeurBloc(b), height: h, rx: 3, fill: "#ece6f5", stroke: encre, "stroke-width": 1 }));
    }
    g.append(el("rect", {
      x: b.x, y: b.y, width: largeurBloc(b), height: h, rx: rond ? 14 : 3,
      fill: sous ? "#f6f2fb" : rond ? "#f7f8fa" : "#fff", stroke: encre, "stroke-width": choisi ? 2.4 : sous ? 2 : 1.6,
    }));
  }
  const lib = libelleBloc(b);
  const taille = lib.length > 22 ? 9 : lib.length > 14 ? 10 : lib.length > 9 ? 12 : 14;
  g.append(el("text", {
    x: b.type === "gain" ? b.x + BLOC_L * 0.36 : b.x + largeurBloc(b) / 2, y: b.y + h / 2 + taille / 3,
    "text-anchor": "middle", "font-size": taille, fill: encre, class: "bloc-texte",
  }, lib));
  g.append(el("text", { x: b.x + largeurBloc(b) / 2, y: b.y + h + (sous ? 17 : 13), "text-anchor": "middle", "font-size": 10, fill: "#98a2b3" }, b.id));
  // ports d'entree (avec le signe d'une somme, ou le nom de l'entree)
  const signes = String(b.signes || "++");
  const noms = nomsEntrees(b);
  for (let k = 0; k < nbEntrees(b); k++) {
    const [x, y] = portEntree(b, k);
    const p = el("circle", { cx: x, cy: y, r: 4.5, fill: "#fff", stroke: "#344054", "stroke-width": 1.2, class: "bloc-port-entree", "data-bloc": b.id, "data-entree": k });
    g.append(p);
    if (b.type === "somme") g.append(el("text", { x: x - 12, y: y - 4, "font-size": 12, fill: encre }, signes[k] || "+"));
    else if (noms[k]) g.append(el("text", { x: x + 7, y: y + 3.5, "font-size": 9, fill: "#667085", class: "bloc-nom-port" }, noms[k]));
  }
  const sorties = nomsSorties(b);
  sorties.forEach((nom, k) => {
    const [x, y] = portSortie(b, k);
    const p = el("circle", { cx: x, cy: y, r: 5, fill: "#344054", class: "bloc-port-sortie", "data-bloc": b.id, "data-sortie": k });
    p.addEventListener("mousedown", (ev) => {
      ev.stopPropagation();
      B.drag = { genre: "lien", de: b.id, sortie: k, fin: null };
    });
    g.append(p);
    if (sorties.length > 1 || sous) g.append(el("text", { x: x - 8, y: y + 3.5, "text-anchor": "end", "font-size": 9, fill: "#667085", class: "bloc-nom-port" }, nom));
  });
  g.addEventListener("mousedown", (ev) => {
    if (ev.target.classList.contains("bloc-port-sortie")) return;
    ev.stopPropagation();
    const s = niveau();
    const [mx, my] = pointSvgBlocs(ev);
    // double-clic sur un sous-systeme : on y entre. Il se detecte sur deux
    // mousedown rapproches, pas sur dblclick : le premier clic redessine le
    // schema, et le bloc vise n'existe plus quand le dblclick arriverait.
    const maintenant = ev.timeStamp;
    if (sous && B.clic && B.clic.id === b.id && maintenant - B.clic.t < 450 && B.chemin.join("/") === B.clic.chemin) {
      B.clic = null;
      B.drag = null;
      entrerSousSysteme(b.id);
      return;
    }
    B.clic = { id: b.id, t: maintenant, chemin: B.chemin.join("/") };
    // Maj+clic : ajoute le bloc a la selection, ou l'en retire
    if (ev.shiftKey) {
      if (B.multi.has(b.id)) B.multi.delete(b.id);
      else B.multi.add(b.id);
      B.sel = B.multi.size ? { genre: "bloc", id: [...B.multi][B.multi.size - 1] } : null;
      B.drag = null;
    } else {
      // un clic sur un bloc deja choisi d'un groupe garde le groupe (pour le deplacer)
      if (!B.multi.has(b.id)) B.multi = new Set([b.id]);
      B.sel = { genre: "bloc", id: b.id };
      const groupe = s.blocs.filter((x) => B.multi.has(x.id));
      B.drag = { genre: "bloc", id: b.id, dx: mx - b.x, dy: my - b.y, bouge: false, groupe: groupe.map((x) => ({ id: x.id, x: x.x, y: x.y })), debut: [mx, my] };
    }
    dessinerSchemaBlocs();
    rendreProprietesBloc();
    majBoutonsGroupe();
  });
  return g;
}

function pointSvgBlocs(ev) {
  const svg = $("#blocs-svg");
  const r = svg.getBoundingClientRect();
  const vb = svg.viewBox.baseVal;
  return [vb.x + ((ev.clientX - r.left) * vb.width) / r.width, vb.y + ((ev.clientY - r.top) * vb.height) / r.height];
}

/** Les liens qui visent un port qu'un sous-systeme n'a plus (port supprime,
 *  renomme...) sont retires, a tous les niveaux. */
function nettoyerLiensSousSystemes(n) {
  const parId = Object.fromEntries(n.blocs.map((b) => [b.id, b]));
  n.liens = n.liens.filter((l) => {
    const a = parId[l.de], c = parId[l.vers];
    if (a && a.type === "sous_systeme" && (+l.sortie || 0) >= nbSorties(a)) return false;
    if (c && c.type === "sous_systeme" && (+l.entree || 0) >= nbEntrees(c)) return false;
    return true;
  });
  for (const b of n.blocs) if (b.type === "sous_systeme" && b.schema) nettoyerLiensSousSystemes(b.schema);
}

/** Le schema change : redessin, et le reseau (qui porte les modeles) est
 *  marque modifie -- c'est ce qui l'enregistre et le rend annulable. */
function schemaModifie() {
  B.compile = null;
  $("#blocs-equations").textContent = "";
  $("#blocs-etat").textContent = "";
  const racine = schemaCourant();
  if (racine) nettoyerLiensSousSystemes(racine);
  const niv = niveau();
  // la selection ne garde que des blocs qui existent encore
  if (niv) B.multi = new Set([...B.multi].filter((id) => niv.blocs.some((b) => b.id === id)));
  rendrePaletteBlocs();
  rendreAriane();
  dessinerSchemaBlocs();
  rendreProprietesBloc();
  rendreParametresModele();
  majBoutonsGroupe();
  modifie();
}

/** Les valeurs d'un nouveau bloc, champ par champ. */
const VALEURS_DEFAUT = {
  K: 1, tau: 0.05, Kp: 1, Ki: 10, Kd: 0, Tf: 0.01, T1: 0.1, T2: 0.05, wn: 100, zeta: 0.7,
  min: -1, max: 1, num: [1], den: [1, 1], valeur: 0, signes: "+-", n: 1, expr: "u",
  Imax: 1.2, Mmax: 1.15, Taw: 0.002, R: 0.01, L: 0.1, C: 0.1,
};
const PREFIXE_ID = {
  premier_ordre: "f", integrateur: "i", derivee: "d", consigne: "c", sortie: "o", entree: "u",
  sous_systeme: "ss", port_entree: "pe", port_sortie: "ps", fonction_transfert: "ft", avance_retard: "ar",
  second_ordre: "so", saturation: "sat", implicite: "imp",
};

function ajouterBloc(type) {
  const s = niveau();
  if (!s) return;
  const spec = specType(type);
  if (!spec) return;
  let k = 1;
  const prefixe = PREFIXE_ID[type] || type[0];
  // l'identifiant est unique dans tout le modele, pas seulement dans ce niveau
  const pris = new Set();
  const collecter = (n) => n.blocs.forEach((x) => { pris.add(x.id); if (x.schema) collecter(x.schema); });
  collecter(schemaCourant());
  while (pris.has(prefixe + k)) k++;
  const b = { id: prefixe + k, type };
  for (const c of spec.champs) {
    if (c === "schema") continue;
    const v = VALEURS_DEFAUT[c] ?? "";
    b[c] = Array.isArray(v) ? [...v] : v;
  }
  if (type === "consigne") b.nom = ["P_ref", "Q_ref"].find((n) => !s.blocs.some((x) => x.nom === n)) || b.id;
  if (type === "entree") b.signal = "V";
  if (type === "sortie") b.signal = ["P", "Q"].find((n) => !s.blocs.some((x) => x.type === "sortie" && x.signal === n)) || "P";
  if (type === "sous_systeme") {
    b.nom = b.id;
    b.schema = { blocs: [], liens: [] };
  }
  if (TYPES_INTERNES.includes(type)) {
    // le port suivant : l'ordre des entrees (ou des sorties) du sous-systeme
    const indices = s.blocs.filter((x) => x.type === type).map((x) => +x.index || 0);
    b.index = indices.length ? Math.max(...indices) + 1 : 0;
    b.nom = (type === "port_entree" ? "u" : "y") + (b.index + 1);
  }
  // une place libre, en balayant la grille
  const occupe = (x, y) => s.blocs.some((o) => Math.abs(o.x - x) < BLOC_L + 10 && y < o.y + hauteurBloc(o) + 20 && o.y < y + 100);
  let x = 40, y = 40;
  while (occupe(x, y)) {
    x += 140;
    if (x > 860) { x = 40; y += 90; }
  }
  Object.assign(b, { x, y });
  s.blocs.push(b);
  B.sel = { genre: "bloc", id: b.id };
  B.multi = new Set([b.id]);
  schemaModifie();
}

function supprimerSelectionBlocs() {
  const s = niveau();
  if (!s || !B.sel) return;
  if (B.sel.genre === "lien") s.liens.splice(B.sel.index, 1);
  else {
    const ids = new Set(B.multi.size ? B.multi : [B.sel.id]);
    s.blocs = s.blocs.filter((b) => !ids.has(b.id));
    s.liens = s.liens.filter((l) => !ids.has(l.de) && !ids.has(l.vers));
  }
  B.sel = null;
  B.multi = new Set();
  schemaModifie();
}

// ----------------------------------------------------------- proprietes
function rendreProprietesBloc() {
  const n = $("#blocs-proprietes");
  n.textContent = "";
  const s = niveau();
  if (!s || !B.sel) {
    n.append(el("p", { class: "doux" }, "Cliquer un bloc ou un lien. Glisser depuis le point noir d'un bloc (sa sortie) jusqu'a l'entree d'un autre pour les relier ; Suppr efface la selection. Maj+clic ou un rectangle tire sur le fond selectionnent plusieurs blocs, a regrouper en sous-systeme."));
    return;
  }
  if (B.sel.genre === "lien") {
    const l = s.liens[B.sel.index];
    if (!l) return;
    const a = s.blocs.find((x) => x.id === l.de);
    const sortie = a && nbSorties(a) > 1 ? `.${nomsSorties(a)[+l.sortie || 0]}` : "";
    n.append(el("p", {}, `lien ${l.de}${sortie} → ${l.vers} (entree ${(+l.entree || 0) + 1})`));
    n.append(el("button", { type: "button", onclick: supprimerSelectionBlocs }, "Supprimer le lien"));
    return;
  }
  if (B.multi.size > 1) {
    n.append(el("h3", {}, `${B.multi.size} blocs selectionnes`));
    n.append(el("p", { class: "doux blocs-aide" }, [...B.multi].join(", ")));
    n.append(el("button", { type: "button", class: "primaire", onclick: regrouperSelection, title: "Les blocs choisis deviennent un seul bloc, avec une entree par signal qui entre et une sortie par signal utilise dehors" }, "Regrouper en sous-systeme"));
    n.append(" ", el("button", { type: "button", onclick: supprimerSelectionBlocs }, "Supprimer les blocs"));
    return;
  }
  const b = s.blocs.find((x) => x.id === B.sel.id);
  const spec = b && specType(b.type);
  if (!b || !spec) return;
  n.append(el("h3", {}, spec.titre));
  const ligne = (lab, champ) => n.append(el("label", { class: "blocs-champ" }, el("span", {}, lab), champ));
  const id = el("input", { value: b.id, size: 8 });
  id.addEventListener("change", () => {
    const nouv = id.value.trim();
    // un identifiant est unique dans tout le modele aplati : ni deja pris, ni de double souligne
    const pris = new Set();
    const collecter = (nv) => nv.blocs.forEach((x) => { pris.add(x.id); if (x.schema) collecter(x.schema); });
    collecter(schemaCourant());
    if (!/^[A-Za-z][A-Za-z0-9_]*$/.test(nouv) || nouv.includes("__") || pris.has(nouv)) {
      etat("identifiant invalide ou deja pris (lettres, chiffres et _ ; pas de double souligne)", true);
      id.value = b.id;
      return;
    }
    for (const l of s.liens) {
      if (l.de === b.id) l.de = nouv;
      if (l.vers === b.id) l.vers = nouv;
    }
    b.id = nouv;
    B.sel.id = nouv;
    B.multi = new Set([nouv]);
    schemaModifie();
  });
  ligne("identifiant", id);
  for (const c of spec.champs) {
    if (c === "schema") continue;
    let champ;
    if (c === "signal") {
      champ = el("select", {});
      for (const v of b.type === "entree" ? B.types.entrees : B.types.sorties)
        champ.append(el("option", { value: v, selected: b[c] === v ? "" : null }, v));
    } else if (LISTES.includes(c)) {
      champ = el("input", { value: (Array.isArray(b[c]) ? b[c] : []).join(", "), size: 12, "data-champ": c });
    } else champ = el("input", { value: b[c] ?? "", size: 12, "data-champ": c });
    champ.addEventListener("change", () => {
      const v = champ.value.trim();
      if (LISTES.includes(c)) b[c] = v.split(",").map((t) => t.trim()).filter((t) => t !== "").map((t) => (Number.isNaN(+t) ? t : +t));
      else if (c === "n") b[c] = Math.max(1, parseInt(v, 10) || 1);
      else if (c === "index") b[c] = Math.max(0, parseInt(v, 10) || 0);
      else b[c] = v !== "" && !Number.isNaN(+v) && !["signes", "expr", "nom"].includes(c) ? +v : v;
      if (c === "signes" || c === "n") {
        const k = nbEntrees(b);
        s.liens = s.liens.filter((l) => l.vers !== b.id || (+l.entree || 0) < k);
      }
      schemaModifie();
    });
    const aide = {
      K: "nombre ou expression des parametres", tau: "s", signes: "ex. +- : une entree par signe",
      expr: b.type === "implicite" ? "residu en u1..un et y : la sortie y verifie expr = 0" : "en u (ou u1, u2...) et les parametres",
      nom: b.type === "consigne" ? "entree de service, calee au load flow" : b.type === "sous_systeme" ? "nom affiche et dans le fil d'Ariane" : "nom du port, affiche sur le sous-systeme",
      index: "rang du port, a partir de 0 : fixe l'ordre des entrees (ou des sorties) du sous-systeme",
      Imax: "pu du convertisseur ; la consigne (d*, q*) est ramenee dans |i| <= Imax, direction conservee",
      Mmax: "indice de modulation maximal (1,15 : modulation vectorielle) ; 0 desactive la saturation",
      Taw: "s ; retro-calcul : relier aw_d / aw_q a l'entree aw des PI (pi_aw) de la boucle en amont",
      Tf: "s ; constante de temps du filtre de la derivee, > 0",
      T1: "s ; avance", T2: "s ; retard, > 0", wn: "rad/s ; pulsation propre", zeta: "amortissement",
      num: "coefficients de s, puissances decroissantes, separes par des virgules (nombres ou parametres)",
      den: "coefficients de s, puissances decroissantes ; degre >= celui du numerateur",
      min: "borne basse", max: "borne haute",
    }[c];
    ligne(c, champ);
    if (aide) n.append(el("div", { class: "doux blocs-aide" }, aide));
  }
  if (b.type === "sous_systeme") {
    const entrees = portsDe(b, "port_entree"), sorties = portsDe(b, "port_sortie");
    n.append(el("p", { class: "doux blocs-aide" }, `entrees : ${entrees.join(", ") || "aucune"} · sorties : ${sorties.join(", ") || "aucune"}`));
    n.append(el("button", { type: "button", class: "primaire", "data-action": "entrer", title: "ou double-clic sur le bloc", onclick: () => entrerSousSysteme(b.id) }, "Entrer"));
    n.append(" ", el("button", { type: "button", "data-action": "degrouper", onclick: degrouperSelection, title: "Remet les blocs du sous-systeme au niveau courant, reliés comme avant" }, "Degrouper"));
    n.append(" ");
  }
  n.append(el("button", { type: "button", onclick: supprimerSelectionBlocs }, "Supprimer le bloc"));
}

/** Les champs dont la valeur est une liste (separee par des virgules). */
const LISTES = ["num", "den"];

function rendreParametresModele() {
  const n = $("#blocs-parametres");
  n.textContent = "";
  const s = schemaCourant();
  if (!s) return;
  s.parametres = s.parametres || {};
  for (const [k, v] of Object.entries(s.parametres)) {
    const val = el("input", { type: "number", step: "any", value: v, style: "width:7em" });
    val.addEventListener("change", () => {
      if (Number.isNaN(parseFloat(val.value))) return;
      s.parametres[k] = parseFloat(val.value);
      schemaModifie();
    });
    const sup = el("button", { type: "button", title: "Retirer ce parametre" }, "✕");
    sup.addEventListener("click", () => {
      delete s.parametres[k];
      schemaModifie();
    });
    n.append(el("div", { class: "blocs-param" }, el("code", {}, k), val, sup));
  }
  const nom = el("input", { placeholder: "nom", size: 6 });
  const val = el("input", { type: "number", step: "any", placeholder: "valeur", style: "width:7em" });
  const aj = el("button", { type: "button" }, "+ parametre");
  aj.addEventListener("click", () => {
    const k = nom.value.trim();
    if (!/^[A-Za-z][A-Za-z0-9_]*$/.test(k)) return etat("nom de parametre invalide", true);
    s.parametres[k] = parseFloat(val.value) || 0;
    schemaModifie();
  });
  n.append(el("div", { class: "blocs-param" }, nom, val, aj));
}

// ------------------------------------------------------------ compilation
async function compilerBlocs() {
  const s = schemaCourant();
  if (!s) return;
  if (s.genre === "equations") return verifierModeleEquations();
  const sortie = $("#blocs-equations");
  sortie.textContent = "";
  $("#blocs-etat").textContent = "compilation...";
  try {
    const d = await api("/api/blocs/compiler", { schema_bloc: s });
    B.compile = d;
    $("#blocs-etat").textContent = `${d.etats.length} etat(s), ${d.algebriques.length} algebrique(s), consignes ${d.consignes.join(", ")}`;
    const latex = $("#blocs-latex").checked;
    const liste = (titre, xs) => sortie.append(el("div", {}, el("b", {}, titre + " : "), xs.length ? xs.join(", ") : "aucun"));
    liste("Parametres", Object.entries(d.parametres).map(([k, v]) => `${k} = ${v}`));
    liste("Etats", d.etats);
    liste("Algebriques", d.algebriques);
    liste("Consignes (calees au load flow)", d.consignes);
    liste("Entrees du reseau", d.entrees_reseau);
    const pre = el("pre", { class: "mono blocs-eq" });
    for (const e of d.equations) pre.append((latex ? e.latex : e.texte) + "\n");
    sortie.append(pre);
  } catch (err) {
    B.compile = null;
    $("#blocs-etat").textContent = "";
    sortie.append(el("div", { class: "pb erreur" }, err.message));
  }
}

// ------------------------------------------------------------ les modeles
function nomLibre(base) {
  const m = modelesBlocs();
  let nom = base, k = 2;
  while (m[nom]) nom = `${base}_${k++}`;
  return nom;
}

/** Ouvre le modele ``nom`` au premier niveau, sans selection. */
function ouvrirModeleBlocs(nom) {
  B.modele = nom;
  B.chemin = [];
  B.sel = null;
  B.multi = new Set();
  B.clic = null;
  B.compile = null;
}

function nouveauModeleBlocs(exemple) {
  if (!E.tables) return etat("charger ou creer un reseau d'abord", true);
  if (!B.types) return etat("les exemples ne sont pas encore charges", true);
  const s = structuredClone(B.types.exemples[exemple].schema);
  s.nom = nomLibre(s.nom);
  modelesBlocs()[s.nom] = s;
  ouvrirModeleBlocs(s.nom);
  modifie();
  rendreBlocs();
}

function renommerModeleBlocs() {
  const s = schemaCourant();
  if (!s) return;
  const nouv = (prompt("Nouveau nom du modele", s.nom) || "").trim();
  if (!nouv || nouv === s.nom) return;
  if (!/^[A-Za-z][A-Za-z0-9_]*$/.test(nouv) || modelesBlocs()[nouv]) return etat("nom invalide ou deja pris", true);
  const m = modelesBlocs();
  delete m[s.nom];
  for (const r of lignes("T_Blocs")) if (r[2] === s.nom) r[2] = nouv;
  s.nom = nouv;
  m[nouv] = s;
  B.modele = nouv;
  modifie();
  rendreBlocs();
}

function dupliquerModeleBlocs() {
  const s = schemaCourant();
  if (!s) return;
  const c = structuredClone(s);
  c.nom = nomLibre(s.nom);
  modelesBlocs()[c.nom] = c;
  ouvrirModeleBlocs(c.nom);
  modifie();
  rendreBlocs();
}

function supprimerModeleBlocs() {
  const s = schemaCourant();
  if (!s) return;
  const n = lignes("T_Blocs").filter((r) => r[2] === s.nom).length;
  if (n) return etat(`${n} element(s) de T_Blocs emploient ${s.nom} : les retirer d'abord`, true);
  delete modelesBlocs()[s.nom];
  B.modele = null;
  B.sel = null;
  B.multi = new Set();
  B.chemin = [];
  modifie();
  rendreBlocs();
}

// ------------------------------------------------------------ evenements
function installerBlocs() {
  const svg = $("#blocs-svg");
  if (!svg) return;
  svg.addEventListener("mousemove", (ev) => {
    if (!B.drag) return;
    const s = niveau();
    const [mx, my] = pointSvgBlocs(ev);
    if (B.drag.genre === "bloc") {
      // tout le groupe suit le bloc saisi
      const px = Math.round((mx - B.drag.dx) / 10) * 10, py = Math.round((my - B.drag.dy) / 10) * 10;
      const ref = B.drag.groupe.find((g) => g.id === B.drag.id);
      for (const g of B.drag.groupe) {
        const b = s.blocs.find((x) => x.id === g.id);
        if (!b) continue;
        b.x = g.x + (px - ref.x);
        b.y = g.y + (py - ref.y);
      }
      B.drag.bouge = true;
    } else B.drag.fin = [mx, my];
    dessinerSchemaBlocs();
  });
  const fin = (ev) => {
    if (!B.drag) return;
    const d = B.drag;
    B.drag = null;
    const s = niveau();
    if (d.genre === "bloc") {
      if (d.bouge) schemaModifie();
      return;
    }
    if (d.genre === "rect") {
      // rectangle de selection : les blocs qu'il touche ; Maj ajoute a l'existant
      if (d.fin && (Math.abs(d.fin[0] - d.debut[0]) > 4 || Math.abs(d.fin[1] - d.debut[1]) > 4)) {
        const [x0, y0] = d.debut, [x1, y1] = d.fin;
        const xa = Math.min(x0, x1), xb = Math.max(x0, x1), ya = Math.min(y0, y1), yb = Math.max(y0, y1);
        const touches = s.blocs.filter((b) => b.x < xb && b.x + largeurBloc(b) > xa && b.y < yb && b.y + hauteurBloc(b) > ya).map((b) => b.id);
        B.multi = new Set(d.ajout ? [...d.base, ...touches] : touches);
        B.sel = B.multi.size ? { genre: "bloc", id: [...B.multi][B.multi.size - 1] } : null;
      }
      dessinerSchemaBlocs();
      rendreProprietesBloc();
      majBoutonsGroupe();
      return;
    }
    // un lien : lache sur une entree (ou sur un bloc, sa premiere entree libre)
    const cible = document.elementFromPoint(ev.clientX, ev.clientY);
    let vers = null, k = 0;
    if (cible && cible.classList.contains("bloc-port-entree")) {
      vers = cible.dataset.bloc;
      k = +cible.dataset.entree;
    } else {
      const g = cible && cible.closest && cible.closest(".bloc");
      if (g) {
        vers = g.dataset.blocId;
        const b = s.blocs.find((x) => x.id === vers);
        const prises = new Set(s.liens.filter((l) => l.vers === vers).map((l) => +l.entree || 0));
        k = [...Array(nbEntrees(b)).keys()].find((i) => !prises.has(i));
        if (k === undefined) vers = null;
      }
    }
    if (!vers || vers === d.de || !nbEntrees(s.blocs.find((x) => x.id === vers))) {
      dessinerSchemaBlocs();
      return;
    }
    s.liens = s.liens.filter((l) => !(l.vers === vers && (+l.entree || 0) === k));
    s.liens.push(d.sortie ? { de: d.de, vers, entree: k, sortie: d.sortie } : { de: d.de, vers, entree: k });
    B.sel = { genre: "lien", index: s.liens.length - 1 };
    B.multi = new Set();
    schemaModifie();
  };
  svg.addEventListener("mouseup", fin);
  svg.addEventListener("mouseleave", (ev) => B.drag && fin(ev));
  // le fond : un clic deselectionne, un glisser tire un rectangle de selection
  svg.addEventListener("mousedown", (ev) => {
    const [mx, my] = pointSvgBlocs(ev);
    const base = ev.shiftKey ? new Set(B.multi) : new Set();
    if (!ev.shiftKey) {
      B.sel = null;
      B.multi = new Set();
    }
    B.drag = { genre: "rect", debut: [mx, my], fin: null, ajout: ev.shiftKey, base };
    dessinerSchemaBlocs();
    rendreProprietesBloc();
    majBoutonsGroupe();
  });
  document.addEventListener("keydown", (ev) => {
    if (!$("#onglet-blocs").classList.contains("actif") || !B.sel) return;
    if (["INPUT", "SELECT", "TEXTAREA"].includes(document.activeElement.tagName)) return;
    if (ev.key === "Delete" || ev.key === "Backspace") {
      ev.preventDefault();
      supprimerSelectionBlocs();
    }
  });
  // Echap remonte d'un niveau (hors saisie)
  document.addEventListener("keydown", (ev) => {
    if (ev.key !== "Escape" || !B.chemin.length || !$("#onglet-blocs").classList.contains("actif")) return;
    if (["INPUT", "SELECT", "TEXTAREA"].includes(document.activeElement.tagName)) return;
    allerAuNiveau(B.chemin.length - 1);
  });
  $("#blocs-modele").addEventListener("change", (ev) => {
    ouvrirModeleBlocs(ev.target.value);
    $("#blocs-equations").textContent = "";
    rendreBlocs();
  });
  const ex = $("#blocs-exemple");
  chargerTypesBlocs().then(() => {
    for (const [k, v] of Object.entries(B.types.exemples)) ex.append(el("option", { value: k }, v.titre));
  }).catch(() => {});
  $("#btn-blocs-nouveau").addEventListener("click", () => nouveauModeleBlocs(ex.value));
  $("#btn-blocs-renommer").addEventListener("click", renommerModeleBlocs);
  $("#btn-blocs-dupliquer").addEventListener("click", dupliquerModeleBlocs);
  $("#btn-blocs-supprimer").addEventListener("click", supprimerModeleBlocs);
  $("#btn-blocs-compiler").addEventListener("click", compilerBlocs);
  $("#btn-blocs-regrouper").addEventListener("click", regrouperSelection);
  $("#btn-blocs-degrouper").addEventListener("click", degrouperSelection);
  $("#btn-blocs-enregistrer").addEventListener("click", enregistrerBibliotheque);
  $("#btn-blocs-biblio").addEventListener("click", basculerBibliotheque);
  $("#blocs-latex").addEventListener("change", () => B.compile && compilerBlocs());
  $("#btn-blocs-poser").addEventListener("click", () => {
    if (!B.modele) return;
    activerOnglet("constructeur");
    outil("Blocs");
    etat(`cliquer le bus ou raccorder un element ${B.modele}`);
  });
}

// --------------------------------------------- appels aux routes de FR-051
/** Un appel (GET, POST ou DELETE) aux routes /api/blocs : comme ``api()``,
 *  mais avec le statut en main, pour dire clairement ce qu'il faut corriger
 *  quand la route manque (404) ou refuse le schema (422). */
async function appelBlocs(methode, chemin, corps) {
  const opt = { method: methode };
  if (corps !== undefined) {
    opt.headers = { "Content-Type": "application/json" };
    opt.body = JSON.stringify(corps);
  }
  let r;
  try {
    r = await fetch(chemin, opt);
  } catch {
    throw new Error("le serveur ne repond pas : verifier qu'il tourne toujours");
  }
  const txt = await r.text();
  let data = null;
  try {
    data = JSON.parse(txt);
  } catch {
    // la reponse n'est pas du JSON : le statut parlera
  }
  if (r.ok) return data;
  const detail = data && data.detail;
  const msg = typeof detail === "string"
    ? detail
    : Array.isArray(detail)
      ? detail.map((d) => `${(d.loc || []).slice(1).join(".") || "requete"} : ${d.msg}`).join(" ; ")
      : txt.slice(0, 200);
  if (r.status === 404 && (!detail || detail === "Not Found"))
    throw new Error("cette fonction n'existe pas dans la version du serveur lancee (route " + chemin.split("?")[0] + " absente) : relancer l'atelier avec la version a jour");
  if (r.status === 422) throw new Error("le serveur refuse ce schema : " + msg);
  throw new Error(msg || `erreur ${r.status}`);
}

/** Les blocs et liens du niveau remplaces par ceux d'un schema que le
 *  serveur a rendu (le niveau garde son objet : le reste du schema, ses
 *  parametres et ses autres niveaux, ne bouge pas). */
function remplacerNiveau(niv, schema) {
  niv.blocs = schema.blocs;
  niv.liens = schema.liens;
}

/** Le schema qu'on envoie au serveur pour regrouper ou degrouper : le modele
 *  lui-meme au premier niveau ; a l'interieur d'un sous-systeme, le contenu
 *  du niveau, habille du nom et des parametres du modele. */
function schemaPourServeur() {
  const racine = schemaCourant(), niv = niveau();
  if (niv === racine) return racine;
  return { nom: racine.nom, parametres: racine.parametres || {}, blocs: niv.blocs, liens: niv.liens };
}

async function regrouperSelection() {
  const niv = niveau();
  if (!niv || B.multi.size < 1) return etat("selectionner d'abord les blocs a regrouper (Maj+clic ou rectangle)", true);
  const nom = (prompt("Nom du sous-systeme", "sous_systeme") || "").trim();
  if (!nom) return;
  const ids = niv.blocs.filter((b) => B.multi.has(b.id)).map((b) => b.id);
  const avant = new Set(niv.blocs.map((b) => b.id));
  const modele = B.modele, chemin = B.chemin.join("/");
  try {
    const d = await appelBlocs("POST", "/api/blocs/regrouper", { schema_bloc: schemaPourServeur(), ids, nom });
    if (B.modele !== modele || B.chemin.join("/") !== chemin) return;
    remplacerNiveau(niv, d.schema);
    const nouveau = niv.blocs.find((b) => !avant.has(b.id) && b.type === "sous_systeme");
    B.multi = new Set(nouveau ? [nouveau.id] : []);
    B.sel = nouveau ? { genre: "bloc", id: nouveau.id } : null;
    schemaModifie();
    etat(`${ids.length} bloc(s) regroupes dans le sous-systeme ${nom}`);
  } catch (err) {
    etat("Regrouper : " + err.message, true);
    $("#blocs-etat").textContent = "Regrouper : " + err.message;
  }
}

async function degrouperSelection() {
  const niv = niveau();
  const b = niv && B.sel && B.sel.genre === "bloc" && niv.blocs.find((x) => x.id === B.sel.id);
  if (!b || b.type !== "sous_systeme") return etat("selectionner un sous-systeme a degrouper", true);
  const modele = B.modele, chemin = B.chemin.join("/");
  const avant = new Set(niv.blocs.map((x) => x.id));
  try {
    const d = await appelBlocs("POST", "/api/blocs/degrouper", { schema_bloc: schemaPourServeur(), id: b.id });
    if (B.modele !== modele || B.chemin.join("/") !== chemin) return;
    remplacerNiveau(niv, d.schema);
    // les blocs qui apparaissent sont ceux du sous-systeme, rendus au niveau
    B.multi = new Set(niv.blocs.filter((x) => !avant.has(x.id)).map((x) => x.id));
    B.sel = B.multi.size ? { genre: "bloc", id: [...B.multi][0] } : null;
    schemaModifie();
    etat(`sous-systeme ${b.nom || b.id} degroupe`);
  } catch (err) {
    etat("Degrouper : " + err.message, true);
    $("#blocs-etat").textContent = "Degrouper : " + err.message;
  }
}

/** Les boutons de la barre : actifs selon la selection. */
function majBoutonsGroupe() {
  const r = $("#btn-blocs-regrouper"), d = $("#btn-blocs-degrouper");
  if (!r || !d) return;
  const niv = schemaCourant() && niveau();
  r.disabled = !niv || B.multi.size < 1;
  const un = niv && B.sel && B.sel.genre === "bloc" && B.multi.size === 1 && niv.blocs.find((x) => x.id === B.sel.id);
  d.disabled = !(un && un.type === "sous_systeme");
}

// ----------------------------------------------- un element en blocs (U4)
/** Les tables dont un element se decrit en blocs (port alternatif). */
const TABLES_EN_BLOCS = /^T_(SMs|GFM\ds|GFL\ds|IBs|Loads)$/;

/** Le schema en blocs de l'element ``table[index]`` du reseau affiche, avec
 *  SES parametres, range comme nouveau modele ``<classe>_<n>`` et ouvert dans
 *  l'onglet Blocs. Avec ``remplacer``, le serveur met ce modele a la place de
 *  l'element (``T_Blocs``, meme bus, meme puissance) : le schema agit alors
 *  sur le reseau, et « Revenir au modele natif » rend l'element.
 *  Rend le nom du modele ; une erreur du serveur est lancee telle quelle. */
async function editerEnBlocs(table, index, remplacer = false) {
  if (!E.tables) throw new Error("charger ou creer un reseau d'abord");
  if (!lignes(table)[index]) throw new Error(`${table} n'a pas de ligne ${index}`);
  if (remplacer) {
    // le SERVEUR fait tout le remplacement : puissance d'une charge reportee
    // a son bus, Sn, origine gardee dans le modele pour revenir au natif
    const d = await appelBlocs("POST", "/api/blocs/remplacer", { tables: E.tables, table, index });
    if (!d || !d.tables || !d.modele) throw new Error("reponse du serveur sans tables");
    appliquerTablesEquations(d.tables);
    if (typeof d.table === "string" && Number.isInteger(d.index))
      suivreReference({ table, index }, { table: d.table, index: d.index });
    else E.sel = null;
    ouvrirModeleBlocs(d.modele);
    modifie();
    return d.modele;
  }
  const d = await appelBlocs("POST", "/api/blocs/depuis_element", { tables: E.tables, table, index });
  const schema = d.schema;
  const base = String((schema.amorce && schema.amorce.classe) || schema.nom || table).replace(/[^A-Za-z0-9_]/g, "");
  const m = modelesBlocs();
  let k = 1;
  while (m[`${base}_${k}`]) k++;
  const nom = `${base}_${k}`;
  schema.nom = nom;
  m[nom] = schema;
  ouvrirModeleBlocs(nom);
  modifie();
  return nom;
}

/** Les deux gestes, depuis un bouton (fiche ou proprietes) : ``sortie`` est
 *  le noeud ou ecrire une erreur, en plus de la ligne d'etat. */
async function boutonEditerEnBlocs(table, index, remplacer, sortie) {
  if (sortie) sortie.textContent = "";
  if (remplacer && typeof confirmerFonctionnementLibre === "function" &&
      !(await confirmerFonctionnementLibre(`${table}[${index}] sera remplace par son modele a blocs.`)))
    return;
  try {
    const nom = await editerEnBlocs(table, index, remplacer);
    if (typeof fermerFiche === "function") fermerFiche();
    activerOnglet("blocs");
    etat(remplacer
      ? `element remplace par son modele a blocs ${nom} (annulable)`
      : `modele ${nom} cree depuis l'element : onglet Blocs`);
  } catch (err) {
    const msg = "Edition en blocs impossible : " + err.message;
    etat(msg, true);
    if (sortie) sortie.append(el("div", { class: "pb erreur" }, msg));
  }
}

/** Les deux boutons « Editer en blocs » et « Remplacer par son modele a
 *  blocs », pour la fiche et les proprietes d'un element. */
function boutonsEnBlocs(table, index, sortie) {
  const edit = el("button", {
    type: "button", class: "btn-en-blocs primaire", "data-action": "editer-en-blocs",
    title: "Ouvre la description en blocs de cet element (avec ses parametres) comme nouveau modele de l'onglet Blocs ; l'element du reseau n'est pas touche",
  }, "Editer en blocs");
  edit.addEventListener("click", () => boutonEditerEnBlocs(table, index, false, sortie));
  const rempl = el("button", {
    type: "button", class: "btn-en-blocs", "data-action": "remplacer-par-blocs",
    title: "Edite en blocs ET met ce modele a la place de l'element dans le reseau : les modifications du schema agissent alors sur les calculs (annulable)",
  }, "Remplacer par son modele a blocs");
  rempl.addEventListener("click", () => boutonEditerEnBlocs(table, index, true, sortie));
  return [edit, " ", rempl];
}

// ------------------------------------------------- bibliotheque (U5)
async function enregistrerBibliotheque() {
  const s = schemaCourant();
  if (!s) return etat("aucun modele ouvert", true);
  try {
    const d = await appelBlocs("POST", "/api/blocs/bibliotheque", { schema_bloc: s });
    etat(`modele enregistre dans la bibliotheque sous le nom ${d.nom}`);
    $("#blocs-etat").textContent = `enregistre dans la bibliotheque : ${d.nom}`;
    if (!$("#blocs-biblio").hidden) await rendreBibliotheque();
  } catch (err) {
    etat("Bibliotheque : " + err.message, true);
    $("#blocs-etat").textContent = "Bibliotheque : " + err.message;
  }
}

function basculerBibliotheque() {
  const p = $("#blocs-biblio");
  p.hidden = !p.hidden;
  if (!p.hidden) rendreBibliotheque();
}

/** La liste des modeles de la bibliotheque, avec Ouvrir et Supprimer. */
async function rendreBibliotheque() {
  const p = $("#blocs-biblio");
  p.textContent = "";
  p.append(el("b", {}, "Bibliotheque de modeles "), el("span", { class: "doux" }, "(partagee entre les reseaux)"));
  let liste;
  try {
    liste = await appelBlocs("GET", "/api/blocs/bibliotheque");
  } catch (err) {
    p.append(el("div", { class: "pb erreur" }, err.message));
    return;
  }
  if (!liste.length) p.append(el("p", { class: "doux" }, "Vide : « Enregistrer dans la bibliotheque » y range le modele ouvert."));
  for (const m of liste) {
    const ouvrir = el("button", { type: "button", "data-ouvrir": m.nom, title: "Ajoute ce modele au reseau et l'ouvre" }, "Ouvrir");
    ouvrir.addEventListener("click", () => ouvrirDepuisBibliotheque(m.nom));
    const sup = el("button", { type: "button", "data-supprimer": m.nom, title: "Retire ce modele de la bibliotheque (les reseaux qui l'ont deja pris le gardent)" }, "✕");
    sup.addEventListener("click", async () => {
      try {
        await appelBlocs("DELETE", "/api/blocs/bibliotheque/" + encodeURIComponent(m.nom));
        await rendreBibliotheque();
      } catch (err) {
        etat("Bibliotheque : " + err.message, true);
      }
    });
    p.append(el("div", { class: "blocs-param", "data-modele-biblio": m.nom },
      el("code", {}, m.titre && m.titre !== m.nom ? `${m.nom} (${m.titre})` : m.nom),
      el("span", { class: "doux" }, m.modifie ? String(m.modifie).slice(0, 16).replace("T", " ") : ""), ouvrir, sup));
  }
}

async function ouvrirDepuisBibliotheque(nom) {
  if (!E.tables) return etat("charger ou creer un reseau d'abord", true);
  try {
    const d = await appelBlocs("GET", "/api/blocs/bibliotheque/" + encodeURIComponent(nom));
    const schema = d.schema;
    schema.nom = nomLibre(schema.nom || nom);
    modelesBlocs()[schema.nom] = schema;
    ouvrirModeleBlocs(schema.nom);
    modifie();
    await rendreBlocs();
    etat(`modele ${schema.nom} ouvert depuis la bibliotheque`);
  } catch (err) {
    etat("Bibliotheque : " + err.message, true);
    $("#blocs-etat").textContent = "Bibliotheque : " + err.message;
  }
}

installerBlocs();
