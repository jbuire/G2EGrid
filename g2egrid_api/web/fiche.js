/* g2egrid -- la FICHE d'un equipement : un double-clic sur un element du
   schema ouvre une page qui decrit sa COMMANDE en schemas-blocs
   d'automatique : sommateurs, blocs de transfert (gain, PI, premier ordre,
   integrateur), saturation, boucles de retour.

   Un schema par boucle, dans l'ordre de la cascade : synchronisation,
   boucles externes, boucle de courant, partie puissance (le procede vu par
   la commande), bus continu. Les schemas decrivent le fonctionnement libre
   (limitation inactive, convertisseur en service) ; l'axe q se deduit de
   l'axe d par symetrie, avec ses propres termes de couplage.

   Sous les schemas, les equations du modele, repliees : celles que
   g2egrid.fiches ecrit et verifie, terme a terme, contre le modele compile.
   Le choix des schemas (``fiche.schema``) vient du serveur, qui connait la
   classe reelle de l'element. */

"use strict";

// --------------------------------------------------------------- ouverture
async function ouvrirFiche(table, index) {
  const fond = $("#fiche");
  const corps = $("#fiche-corps");
  corps.replaceChildren(el("p", { class: "attente" }, "calcul de la fiche..."));
  $("#fiche-titre").textContent = `${table} [${index}]`;
  fond.hidden = false;
  document.body.classList.add("fiche-ouverte");
  try {
    // les tables de la VUE courante : c'est la que vit la ligne visee
    const d = await api("/api/fiche", { tables: E.tables, table, index });
    E.fiche = d;
    E.ficheRef = { table, index };
    rendreFiche(d);
  } catch (err) {
    corps.replaceChildren(el("div", { class: "pb erreur" }, err.message));
  }
}

function fermerFiche() {
  $("#fiche").hidden = true;
  document.body.classList.remove("fiche-ouverte");
}

function initFiche() {
  $("#fiche-fermer").addEventListener("click", fermerFiche);
  $("#fiche").addEventListener("click", (ev) => {
    if (ev.target.id === "fiche") fermerFiche();
  });
  document.addEventListener("keydown", (ev) => {
    if (ev.key === "Escape" && !$("#fiche").hidden) fermerFiche();
  });
}

// ------------------------------------------------------------------- rendu
function rendreFiche(d) {
  const e = d.element;
  $("#fiche-titre").textContent = `${e.nom} — ${e.modele}`;
  const corps = $("#fiche-corps");
  corps.replaceChildren();
  const info = [`bus ${e.bus.join(" / ")}`];
  if (e.Sn) info.push(`Sn = ${fmt(e.Sn / 1e6, 1)} MVA`);
  info.push(`${d.etats.length} etat(s)`);
  corps.append(el("p", { class: "doux" }, info.join(" · ")));
  if (d.resume) corps.append(el("p", {}, d.resume));
  // U4 : cet element, en blocs editables
  const ref = E.ficheRef;
  if (ref && typeof boutonsEnBlocs === "function") {
    const sortie = el("div", {});
    corps.append(el("div", { class: "fiche-en-blocs" },
      ...boutonsEnBlocs(ref.table, ref.index, sortie),
      el("span", { class: "doux" }, " decrit cet element en blocs, avec ses parametres, dans l'onglet Blocs")), sortie);
  }

  if (ref && typeof sectionEquations === "function" && TABLES_EQUATIONS.test(ref.table))
    corps.append(sectionEquations(ref.table, ref.index, "fiche"));

  const schemas = d.schema && SCHEMAS_BLOCS[d.schema];
  if (schemas) {
    corps.append(el("p", { class: "doux fiche-hyp" },
      "Schemas-blocs de la commande, en fonctionnement libre (limitation de courant inactive, " +
      "convertisseur en service). Grandeurs en pu, repere dq du convertisseur ; l'axe q suit " +
      "l'axe d, avec ses propres termes de couplage."));
    schemas.forEach((sch, i) => {
      const sec = el("section", { class: "fiche-bloc", "data-schema-bloc": sch.titre });
      sec.append(el("h4", {}, `${i + 1}. ${sch.titre}`));
      if (sch.texte) sec.append(el("p", { class: "doux" }, sch.texte));
      sec.append(el("div", { class: "fiche-schema" }, schemaBloc(sch)));
      corps.append(sec);
    });
  } else {
    corps.append(el("p", { class: "doux fiche-hyp" },
      "Pas de schema-bloc pour ce modele : ses equations seules, ci-dessous."));
  }

  corps.append(blocEquationsFiche(d, !schemas));
  corps.append(tableSymboles("Symboles", d.legende, (s) => [s.genre, s.description]));
  if (d.parametres && d.parametres.length)
    corps.append(tableSymboles("Valeurs des parametres", d.parametres, (s) => [fmt(s.valeur, 6), s.description]));
}

/** Les equations du modele, par sous-systeme, repliees sous les schemas. */
function blocEquationsFiche(d, ouvert) {
  const det = el("details", { class: "avance fiche-equations", open: ouvert ? "" : null },
    el("summary", {}, "Equations du modele (verifiees contre le modele compile)"));
  if (d.hypotheses) det.append(el("p", { class: "doux fiche-hyp" }, d.hypotheses));
  const genres = { def: "signal", alg: "algebrique", etat: "etat", sortie: "sortie", residu: "residu" };
  d.blocs.forEach((b, i) => {
    const sec = el("section", { class: "fiche-bloc" });
    sec.append(el("h4", {}, `${i + 1}. ${b.titre}`));
    if (b.texte) sec.append(el("p", { class: "doux" }, b.texte));
    for (const l of b.lignes) {
      const eq = el("div", { class: `fiche-eq fiche-${l.genre}`, title: l.texte });
      const m = el("div", { class: "fiche-math" });
      m.innerHTML = l.mathml; // MathML ecrit par SymPy cote serveur
      const etiq = [genres[l.genre] || l.genre];
      if (l.genre === "etat" && d.equilibre && l.nom in d.equilibre)
        etiq.push(`x0 = ${fmt(d.equilibre[l.nom], 5)}`);
      eq.append(m, el("span", { class: "fiche-genre" }, etiq.join(" · ")));
      sec.append(eq);
    }
    det.append(sec);
  });
  const latex = d.blocs.flatMap((b) => [`% ${b.titre}`, ...b.lignes.map((l) => l.latex)]).join("\n");
  const copier = el("button", {}, "Copier en LaTeX");
  copier.addEventListener("click", () => navigator.clipboard && navigator.clipboard.writeText(latex));
  det.append(copier);
  return det;
}

function tableSymboles(titre, lignes, colonnes) {
  const det = el("details", { class: "avance" }, el("summary", {}, titre));
  const t = el("table", { class: "fiche-table" });
  for (const s of lignes) {
    const td = el("td", {});
    td.innerHTML = s.mathml;
    t.append(el("tr", {}, td, el("td", { class: "mono" }, s.nom), ...colonnes(s).map((c) => el("td", {}, c))));
  }
  det.append(t);
  return det;
}

// ------------------------------------------------------------ textes SVG
/** Un texte SVG dont les ``x_y`` deviennent de vrais indices (tspan
 *  decales, que tous les navigateurs savent placer). */
function texteIndice(attrs, t) {
  const n = el("text", attrs);
  const re = /([^\s_()·,/+−=*|²√]+)_([A-Za-z0-9]+(?:,[A-Za-z]+)?)/g;
  let i = 0;
  let m;
  while ((m = re.exec(t))) {
    if (m.index > i) n.append(document.createTextNode(t.slice(i, m.index)));
    n.append(document.createTextNode(m[1]));
    const sub = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
    sub.setAttribute("dy", "3");
    sub.setAttribute("font-size", "75%");
    sub.textContent = m[2];
    const ret = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
    ret.setAttribute("dy", "-3");
    ret.textContent = "​";
    n.append(sub, ret);
    i = m.index + m[0].length;
  }
  if (i < t.length) n.append(document.createTextNode(t.slice(i)));
  return n;
}

/** Largeur approchee d'un texte (les indices comptent moins). */
const largeurTexte = (t, px = 7.2) =>
  t.replace(/_([A-Za-z0-9]+)/g, (_, s) => "x".repeat(Math.ceil(s.length * 0.7))).length * px;

// ----------------------------------------------------- schemas-blocs
/* Un schema est une CHAINE, de gauche a droite, et des RETOURS.

   Elements de la chaine :
     { t: "e", lbl }                  entree (texte, fleche sortante)
     { t: "o", lbl }                  sortie (fleche entrante, texte)
     { t: "s", g, haut, bas }         sommateur ; g : signe de l'entree de
                                      gauche (+ par defaut) ; haut / bas :
                                      { s: "+" | "−", lbl } entrees laterales
     { t: "b", num, den }             transfert en fraction num / den
     { t: "b", txt }                  bloc a texte (gain, operation)
     { t: "sat", lbl }                saturation
   Chaque element peut porter ``sig`` : le nom du signal qui en sort.

   Retours : { de, vers, lbl, s, n } -- de la sortie de l'element ``de`` au
   bas de l'element ``vers`` (un sommateur, avec le signe ``s``, − par
   defaut), sur le niveau ``n`` (1, 2...) sous la chaine. */
const Y0 = 70; // ligne de la chaine
const R_SOMME = 11;
const ECART = 54; // fil entre deux elements

function largeurElement(c) {
  switch (c.t) {
    case "e":
    case "o":
      return Math.max(28, largeurTexte(c.lbl)) + 6;
    case "s":
      return 2 * R_SOMME;
    case "sat":
      return 54;
    default: {
      const t = c.txt ? [c.txt] : [c.num, c.den];
      return Math.max(46, ...t.map((x) => largeurTexte(x, 7.6))) + 22;
    }
  }
}

function schemaBloc(sch) {
  const chaine = sch.chaine;
  const retours = sch.retours || [];
  const aBas = chaine.some((c) => c.bas);
  const nmax = Math.max(0, ...retours.map((r) => r.n || 1));
  const pos = [];
  let x = 10;
  chaine.forEach((c, i) => {
    const w = largeurElement(c);
    pos.push({ x, w, cx: x + w / 2 });
    x += w + (i < chaine.length - 1 ? ECART : 0);
  });
  const W = x + 10;
  const H = Y0 + (aBas ? 62 : 34) + 26 * nmax;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, class: "fiche-svg schema-bloc", role: "img" });
  const T = (attrs, t) => svg.append(texteIndice(attrs, t));
  const ligne = (d, cls = "lien") => svg.append(el("path", { d, class: cls }));
  const pointe = (x2, y2, dir) => {
    const p = { d: [[-7, -4], [-7, 4]], h: [[-4, 7], [4, 7]], b: [[-4, -7], [4, -7]] }[dir];
    ligne(`M ${x2 + p[0][0]} ${y2 + p[0][1]} L ${x2} ${y2} L ${x2 + p[1][0]} ${y2 + p[1][1]}`);
  };

  // les fils de la chaine, avec le nom du signal
  for (let i = 0; i < chaine.length - 1; i++) {
    const x1 = pos[i].x + pos[i].w;
    const x2 = pos[i + 1].x;
    ligne(`M ${x1} ${Y0} L ${x2} ${Y0}`);
    pointe(x2, Y0, "d");
    if (chaine[i].sig) T({ x: (x1 + x2) / 2, y: Y0 - 7, "text-anchor": "middle", class: "signal" }, chaine[i].sig);
  }

  chaine.forEach((c, i) => {
    const { x, w, cx } = pos[i];
    if (c.t === "e") T({ x: x + w - 4, y: Y0 + 4, "text-anchor": "end", class: "es" }, c.lbl);
    else if (c.t === "o") T({ x: x + 4, y: Y0 + 4, "text-anchor": "start", class: "es" }, c.lbl);
    else if (c.t === "s") {
      svg.append(el("circle", { cx, cy: Y0, r: R_SOMME, class: "somme" }));
      const k = R_SOMME * 0.707;
      ligne(`M ${cx - k} ${Y0 - k} L ${cx + k} ${Y0 + k} M ${cx - k} ${Y0 + k} L ${cx + k} ${Y0 - k}`, "croix");
      T({ x: cx - R_SOMME - 3, y: Y0 - 6, "text-anchor": "end", class: "signe" }, c.g || "+");
      if (c.haut) {
        ligne(`M ${cx} ${Y0 - 40} L ${cx} ${Y0 - R_SOMME}`);
        pointe(cx, Y0 - R_SOMME, "b");
        T({ x: cx, y: Y0 - 45, "text-anchor": "middle", class: "es" }, c.haut.lbl);
        T({ x: cx + 5, y: Y0 - R_SOMME - 4, "text-anchor": "start", class: "signe" }, c.haut.s || "+");
      }
      if (c.bas) {
        ligne(`M ${cx} ${Y0 + 40} L ${cx} ${Y0 + R_SOMME}`);
        pointe(cx, Y0 + R_SOMME, "h");
        T({ x: cx, y: Y0 + 54, "text-anchor": "middle", class: "es" }, c.bas.lbl);
        T({ x: cx + 5, y: Y0 + R_SOMME + 12, "text-anchor": "start", class: "signe" }, c.bas.s || "+");
      }
    } else if (c.t === "sat") {
      svg.append(el("rect", { x, y: Y0 - 20, width: w, height: 40, class: "bloc-tf" }));
      ligne(`M ${x + 8} ${Y0 + 10} L ${x + 18} ${Y0 + 10} L ${x + w - 18} ${Y0 - 10} L ${x + w - 8} ${Y0 - 10}`, "croix");
      if (c.lbl) T({ x: cx, y: Y0 - 25, "text-anchor": "middle", class: "es" }, c.lbl);
    } else {
      const h = c.txt ? 32 : 46;
      svg.append(el("rect", { x, y: Y0 - h / 2, width: w, height: h, class: `bloc-tf ${c.k ? "tf-" + c.k : ""}` }));
      if (c.txt) T({ x: cx, y: Y0 + 4, "text-anchor": "middle", class: "tf" }, c.txt);
      else {
        T({ x: cx, y: Y0 - 6, "text-anchor": "middle", class: "tf" }, c.num);
        ligne(`M ${x + 8} ${Y0} L ${x + w - 8} ${Y0}`, "barre-frac");
        T({ x: cx, y: Y0 + 16, "text-anchor": "middle", class: "tf" }, c.den);
      }
    }
  });

  // les retours, sous la chaine
  for (const r of retours) {
    const n = r.n || 1;
    const yr = Y0 + (aBas ? 48 : 20) + 26 * n;
    const a = pos[r.de];
    const b = pos[r.vers];
    const xa = a.x + a.w + 14 + 8 * (n - 1);
    svg.append(el("circle", { cx: xa, cy: Y0, r: 2.6, class: "noeud" }));
    const bas = chaine[r.vers].t === "s" ? Y0 + R_SOMME : Y0 + (chaine[r.vers].txt ? 16 : 23);
    ligne(`M ${xa} ${Y0} L ${xa} ${yr} L ${b.cx} ${yr} L ${b.cx} ${bas}`);
    pointe(b.cx, bas, "h");
    if (r.lbl) T({ x: xa - 6, y: yr - 5, "text-anchor": "end", class: "signal" }, r.lbl);
    if (chaine[r.vers].t === "s")
      T({ x: b.cx - 6, y: bas + 13, "text-anchor": "end", class: "signe" }, r.s || "−");
  }
  return svg;
}

// ------------------------------------------------- briques des schemas
const E_ = (lbl, sig) => ({ t: "e", lbl, sig });
const O_ = (lbl) => ({ t: "o", lbl });
const S_ = (opts = {}, sig) => ({ t: "s", ...opts, sig });
const PI_ = (kp, ki, sig) => ({ t: "b", num: `${kp} s + ${ki}`, den: "s", sig, k: "pi" });
const TF = (num, den, sig, k) => ({ t: "b", num, den, sig, k });
const B_ = (txt, sig, k) => ({ t: "b", txt, sig, k });
const SAT = (sig) => ({ t: "sat", lbl: "|i*| ≤ I_max", sig });
const INT_WB = (sig) => TF("ω_b", "s", sig, "int");

/** La boucle de courant (axe d) et l'inductance cote onduleur. ``vdc`` :
 *  "1" (bus continu ideal) ou "v_dc" ; ``v`` : la tension vue par la
 *  boucle (v_ed derriere un filtre LCL, v_gd pour un MMC). */
function boucleCourant({ v = "v_ed", w = "ω", vdc = "v_dc", mmc = false } = {}) {
  const chaine = [
    E_("i_sd*"),
    SAT("i_sd*"),
    S_(),
    PI_("K_pCL", "K_iCL"),
    S_({ haut: { s: "+", lbl: `K_ffv·${v}` }, bas: { s: "−", lbl: "ω_ff·L_1·i_sq" } }, "e_d"),
  ];
  if (!mmc && vdc !== "1") chaine.push(B_("÷ v_dc", "m_d"), B_("× v_dc", "e_d", "proc"));
  chaine.push(
    S_({ haut: { s: "−", lbl: v }, bas: { s: "+", lbl: `${w}·L_1·i_sq` } }),
    TF("ω_b", "L_1 s + ω_b R_1", "i_sd", "proc"),
    O_("i_sd"),
  );
  const iProc = chaine.length - 2;
  return {
    titre: "Boucle de courant (axe d)",
    texte:
      "Saturation du module de la consigne, PI, anticipation de la tension et decouplage des axes" +
      (mmc ? " ; les bras fabriquent e_d (modulation compensee)." : vdc === "1" ? "." : " ; la modulation divise par v_dc, l'onduleur la remultiplie.") +
      " Le procede est l'inductance R_1-L_1.",
    chaine,
    retours: [{ de: iProc, vers: 2, lbl: "i_sd" }],
  };
}

const FILTRE_LCL = (w = "ω") => ({
  titre: "Filtre LCL (procede, axe d)",
  texte: "Le condensateur integre la difference des courants ; l'inductance cote reseau porte i_g vers le bus.",
  chaine: [
    E_("i_sd"),
    S_({ haut: { s: "+", lbl: `${w}·C_f·v_eq` } }),
    TF("ω_b", "C_f s", "v_ed", "proc"),
    S_({ haut: { s: "+", lbl: `${w}·L_2·i_gq` }, bas: { s: "−", lbl: "v_gd (bus)" } }),
    TF("ω_b", "L_2 s + ω_b R_2", "i_gd", "proc"),
    O_("i_gd → bus"),
  ],
  retours: [{ de: 4, vers: 1, lbl: "i_gd" }],
});

const PLL = (v = "v_eq") => ({
  titre: "PLL",
  texte: `La PLL annule ${v} : un PI donne la pulsation, un integrateur l'angle, qui referme la boucle par la transformation de Park.`,
  chaine: [
    E_("v_abc"),
    B_("Park (θ_pll)", v, "mes"),
    PI_("K_ppll", "K_ipll"),
    S_({ haut: { s: "+", lbl: "ω_ff" } }, "ω_pll"),
    S_({ haut: { s: "−", lbl: "ω_ref" } }),
    INT_WB("θ_pll"),
    O_("θ_pll"),
  ],
  retours: [{ de: 5, vers: 1, lbl: "θ_pll" }],
});

const EXT_Q = {
  titre: "Boucle externe de puissance reactive",
  texte: "Un PI fixe le courant reactif ; un i_sq negatif fournit du reactif.",
  chaine: [E_("q_e"), S_({ bas: { s: "−", lbl: "q*" } }), PI_("K_pq", "K_iq", "i_sq*"), O_("i_sq*")],
};

const EXT_P = {
  titre: "Boucle externe de puissance active",
  chaine: [E_("p*"), S_({ bas: { s: "−", lbl: "p_e" } }), PI_("K_pp", "K_ip", "i_sd*"), O_("i_sd*")],
};

const EXT_VDC = (texte) => ({
  titre: "Boucle externe de tension continue",
  texte,
  chaine: [E_("v_dc"), S_({ bas: { s: "−", lbl: "v_dc*" } }), PI_("K_pdc", "K_idc", "i_sd*"), O_("i_sd*")],
});

const BUS_DC_GFL = {
  titre: "Bus continu (procede)",
  texte: "La source suit sa consigne au premier ordre ; le condensateur integre la difference entre son courant, celui qu'appelle l'onduleur et la fuite.",
  chaine: [
    E_("i_dc*"),
    TF("1", "1 + T_dc s", "i_dc", "proc"),
    S_({ haut: { s: "−", lbl: "m_d i_sd + m_q i_sq" } }),
    TF("ω_b", "C_dc s", "v_dc", "proc"),
    O_("v_dc"),
  ],
  retours: [{ de: 3, vers: 2, lbl: "G_dc·v_dc" }],
};

const STATISME_PF = (ordre2) => ({
  titre: ordre2 ? "Synchronisation : statisme P-f du 2e ordre" : "Synchronisation : statisme P-f",
  texte: ordre2
    ? "Puissance filtree, statisme m_p filtre a ω_c (inertie virtuelle) ; l'ecart de pulsation, integre, donne l'angle."
    : "Puissance filtree, statisme m_p ; l'ecart de pulsation, integre, donne l'angle du formeur.",
  chaine: [
    E_("p_e"),
    TF("ω_fp", "s + ω_fp", "p_m", "mes"),
    S_({ g: "−", haut: { s: "+", lbl: "p*" } }),
    ordre2 ? TF("m_p ω_c", "s + ω_c", "Δω") : B_("m_p", "Δω"),
    INT_WB("θ"),
    O_("θ   (ω = ω_ref + Δω)"),
  ],
});

const QV_GFM = {
  titre: "Statisme Q-V",
  texte: "La consigne de tension du condensateur baisse quand le reactif depasse sa consigne.",
  chaine: [
    E_("q_e"),
    TF("ω_fq", "s + ω_fq", "q_m", "mes"),
    S_({ g: "−", haut: { s: "+", lbl: "q*" } }),
    B_("n_q"),
    S_({ haut: { s: "+", lbl: "v_e*" } }, "v_ed*"),
    O_("v_ed*"),
  ],
};

const BOUCLE_TENSION = (ref = "v_ed*") => ({
  titre: "Boucle de tension (axe d)",
  texte: "PI sur la tension du condensateur, anticipation du courant reseau et decouplage ; la boucle de courant fermee (schema suivant) y est un bloc.",
  chaine: [
    E_(ref),
    S_(),
    PI_("K_pVL", "K_iVL"),
    S_({ haut: { s: "+", lbl: "K_ffi·i_gd" }, bas: { s: "−", lbl: "ω_ff·C_f·v_eq" } }, "i_sd*"),
    B_("boucle de courant", "i_sd", "cour"),
    S_({ haut: { s: "+", lbl: "ω·C_f·v_eq" }, bas: { s: "−", lbl: "i_gd" } }),
    TF("ω_b", "C_f s", "v_ed", "proc"),
    O_("v_ed"),
  ],
  retours: [{ de: 6, vers: 1, lbl: "v_ed", n: 1 }],
});

const DC_GFM = {
  titre: "Bus continu et sa source",
  texte: "La source continue regle v_dc par un gain proportionnel autour de la puissance programmee, au premier ordre T_dc.",
  chaine: [
    E_("v_dc*"),
    S_(),
    B_("K_pdc"),
    S_({ haut: { s: "+", lbl: "p*/v_dc*" } }, "i_dc*"),
    TF("1", "1 + T_dc s", "i_dc", "proc"),
    S_({ haut: { s: "−", lbl: "m_d i_sd + m_q i_sq" } }),
    TF("ω_b", "C_dc s", "v_dc", "proc"),
    O_("v_dc"),
  ],
  retours: [
    { de: 6, vers: 5, lbl: "G_dc·v_dc", n: 1 },
    { de: 6, vers: 1, lbl: "v_dc", n: 2 },
  ],
};

const MATCHING = {
  titre: "Synchronisation : appariement (matching)",
  texte: "La pulsation est l'image de la tension continue : le condensateur continu tient le role de l'inertie.",
  chaine: [E_("v_dc"), B_("ω_ref / v_dc*", "ω"), S_({ haut: { s: "−", lbl: "ω_ref" } }), INT_WB("θ"), O_("θ")],
};

const DVOC = [
  {
    titre: "Oscillateur virtuel (dVOC) : pulsation",
    texte: "Puissance active normalisee par l'amplitude V_r de l'oscillateur.",
    chaine: [
      E_("p_e"),
      TF("ω_fp", "s + ω_fp", "p_m", "mes"),
      B_("η / V_r²"),
      S_({ g: "−", haut: { s: "+", lbl: "η·p*/v_e*²" } }, "Δω"),
      INT_WB("θ"),
      O_("θ"),
    ],
  },
  {
    titre: "Oscillateur virtuel (dVOC) : amplitude",
    texte: "L'amplitude V_r devient la consigne de la boucle de tension (v_ed* = V_r).",
    chaine: [
      E_("q_e"),
      TF("ω_fq", "s + ω_fq", "q_m", "mes"),
      B_("1 / V_r²"),
      S_({ g: "−", haut: { s: "+", lbl: "q*/v_e*²" }, bas: { s: "+", lbl: "α(v_e*² − V_r²)/v_e*²" } }),
      B_("ω_b·η·V_r"),
      TF("1", "s", "V_r", "int"),
      O_("V_r = v_ed*"),
    ],
  },
];

// --- MMC
const PLL_MMC = { ...PLL("v_gq"), texte: "La PLL annule v_gq, la tension du bus sur l'axe q." };
const EXT_VDC_MMC = {
  titre: "Boucle externe de tension continue",
  texte: "Sur la mesure FILTREE de v_dc ; sans integrateur (K_idc = 0), c'est un statisme de pente K_pdc.",
  chaine: [E_("v_dc"), TF("1", "1 + T_vdc s", "v_dc,f", "mes"), S_({ bas: { s: "−", lbl: "v_dc*" } }), PI_("K_pdc", "K_idc", "i_sd*"), O_("i_sd*")],
};
const ENERGIE_LIBRE = [
  {
    titre: "Energie des sous-modules (libre)",
    texte: "L'energie integre l'ecart des puissances continue et alternative ; la tension des bras en decoule.",
    chaine: [
      E_("v_z·i_dc"),
      S_({ bas: { s: "−", lbl: "p_ac = e_d i_sd + e_q i_sq" } }),
      TF("1", "s", "W", "int"),
      B_("√(W / W_0)", "v_z"),
      O_("v_z"),
    ],
  },
  {
    titre: "Courant continu (procede)",
    texte: "L_dc, R_dc entre le noeud continu et la tension des bras.",
    chaine: [E_("v_dc"), S_({ haut: { s: "−", lbl: "v_z" } }), TF("ω_b", "L_dc s", "i_dc", "proc"), O_("i_dc")],
    retours: [{ de: 2, vers: 1, lbl: "R_dc·i_dc" }],
  },
];
const ENERGIE_REGULEE = [
  {
    titre: "Boucles d'energie et de courant continu",
    texte: "Le PI d'energie, anticipe par la puissance mesuree, fixe i_dc* ; le PI de courant continu fabrique u_z = v_dc − v_z.",
    chaine: [
      E_("W_0"),
      S_(),
      PI_("K_pW", "K_iW"),
      S_({ haut: { s: "+", lbl: "p_e" } }, "i_dc*"),
      S_(),
      PI_("K_pZ", "K_iZ", "u_z"),
      S_(),
      TF("ω_b", "L_dc s", "i_dc", "proc"),
      O_("i_dc"),
    ],
    retours: [
      { de: 7, vers: 6, lbl: "R_dc·i_dc", n: 1 },
      { de: 7, vers: 4, lbl: "i_dc", n: 2 },
    ],
  },
  {
    titre: "Energie des sous-modules",
    chaine: [E_("v_z·i_dc"), S_({ bas: { s: "−", lbl: "p_ac = e_d i_sd + e_q i_sq" } }), TF("1", "s", "W", "int"), O_("W → boucle d'energie")],
  },
];

const SCHEMAS_BLOCS = {
  gfl_pq: [PLL(), EXT_P, EXT_Q, boucleCourant({ w: "ω_pll", vdc: "1" }), FILTRE_LCL("ω_pll")],
  gfl_vdc: [
    PLL(),
    EXT_VDC("La boucle de tension continue fixe le courant actif : elle renvoie au reseau ce que le bus continu recoit."),
    EXT_Q,
    boucleCourant({ w: "ω_pll" }),
    FILTRE_LCL("ω_pll"),
    BUS_DC_GFL,
  ],
  gfl_vdc_port: [
    PLL(),
    EXT_VDC("Sans integrateur (K_idc = 0), c'est un statisme de pente K_pdc ; v_dc est la tension du noeud continu."),
    EXT_Q,
    boucleCourant({ w: "ω_pll" }),
    FILTRE_LCL("ω_pll"),
  ],
  gfm: [STATISME_PF(false), QV_GFM, BOUCLE_TENSION(), boucleCourant({ vdc: "1" }), FILTRE_LCL()],
  gfm2: [STATISME_PF(true), QV_GFM, BOUCLE_TENSION(), boucleCourant({ vdc: "1" }), FILTRE_LCL()],
  gfm_dc: [STATISME_PF(false), QV_GFM, BOUCLE_TENSION(), boucleCourant(), FILTRE_LCL(), DC_GFM],
  gfm2_dc: [STATISME_PF(true), QV_GFM, BOUCLE_TENSION(), boucleCourant(), FILTRE_LCL(), DC_GFM],
  gfm_port: [STATISME_PF(false), QV_GFM, BOUCLE_TENSION(), boucleCourant(), FILTRE_LCL()],
  gfm2_port: [STATISME_PF(true), QV_GFM, BOUCLE_TENSION(), boucleCourant(), FILTRE_LCL()],
  gfm_dvoc: [...DVOC, BOUCLE_TENSION("V_r"), boucleCourant({ vdc: "1" }), FILTRE_LCL()],
  gfm_matching: [MATCHING, QV_GFM, BOUCLE_TENSION(), boucleCourant(), FILTRE_LCL(), DC_GFM],
  mmc: [PLL_MMC, EXT_VDC_MMC, EXT_Q, boucleCourant({ v: "v_gd", w: "ω_pll", mmc: true }), ...ENERGIE_LIBRE],
  mmc_energie: [PLL_MMC, EXT_VDC_MMC, EXT_Q, boucleCourant({ v: "v_gd", w: "ω_pll", mmc: true }), ...ENERGIE_REGULEE],
};

initFiche();
