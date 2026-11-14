/* g2egrid -- l'interface en deux langues (FR / EN).

   L'atelier est ECRIT en francais : ses textes sont dans index.html et dans
   les scripts, la ou ils servent. Ce fichier n'y change rien. Quand la langue
   choisie est l'anglais, il traduit ce que la page AFFICHE :

   * les noeuds de texte, a mesure qu'ils entrent dans la page (un
     MutationObserver), et les attributs title, placeholder et aria-label ;
   * les boites alert / confirm / prompt.

   Le dictionnaire (i18n_en.js) va du texte francais au texte anglais. Trois
   sortes d'entrees, essayees dans cet ordre :

   1. le texte exact ;
   2. un GABARIT a trous, "bus {0} a {1} s" : les trous prennent n'importe
      quoi (un nombre, un nom), et leur contenu est traduit a son tour s'il
      est connu ;
   3. a defaut, le texte est decoupe sur ses separateurs forts (retour a la
      ligne, " -- ", " ; ", " : ") et chaque morceau est traduit s'il est
      connu -- « Validation : <message du serveur> » garde son message.

   Un texte inconnu reste en francais : rien ne casse, il manque une entree.
   I18N.manquants() liste ceux que la page a montres (pour completer le
   dictionnaire). Le texte francais d'origine est garde sur chaque noeud :
   revenir au francais le rend tel quel. Les donnees de l'utilisateur (noms,
   valeurs des champs) ne sont jamais des noeuds de texte de la page ou ne
   sont pas dans le dictionnaire : elles ne bougent pas. */

"use strict";

const I18N = (() => {
  const CLE = "g2egrid.langue";
  const LANGUES = ["fr", "en"];
  const ATTRIBUTS = ["title", "placeholder", "aria-label"];
  // sous ces balises, rien n'est du texte d'interface
  const HORS = new Set(["SCRIPT", "STYLE", "TEXTAREA", "MATH", "CODE", "PRE"]);
  const LETTRES = /[A-Za-zÀ-ÿ]{2}/;
  const SEPARATEURS = ["\n", " -- ", " · ", " ; ", " : "];
  const LONGUEUR_MAX = 600;

  let langue = "fr";
  const exact = new Map(); // texte francais normalise -> anglais
  let gabarits = []; // [{re, ancre, en, poids}]
  const cache = new Map(); // texte normalise -> anglais | null
  const inconnus = new Set();
  let observateur = null;

  /** Le texte tel que le dictionnaire le range : sans espaces de bord, les
   *  retours a la ligne debarrasses de l'indentation du source. */
  const normaliser = (s) => s.replace(/[ \t\r]*\n[ \t\r]*/g, "\n").trim();

  const echapper = (s) => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");

  /** Charge un dictionnaire {francais: anglais}. Une cle qui porte des trous
   *  {0}, {1}... devient un gabarit. */
  function charger(dico) {
    for (const [fr, en] of Object.entries(dico)) {
      const cle = normaliser(fr);
      if (!/\{\d+\}/.test(cle)) {
        exact.set(cle, en);
        continue;
      }
      const morceaux = cle.split(/(\{\d+\})/);
      const ordre = [];
      let source = "", ancre = "", poids = 0;
      for (const m of morceaux) {
        const trou = /^\{(\d+)\}$/.exec(m);
        if (trou) {
          ordre.push(+trou[1]);
          source += "([\\s\\S]*?)";
        } else {
          source += echapper(m);
          poids += m.length;
          if (m.length > ancre.length) ancre = m;
        }
      }
      if (!LETTRES.test(ancre)) continue; // "{0} ({1})" : rien a traduire
      gabarits.push({ re: new RegExp("^" + source + "$"), ancre, en, ordre, poids });
    }
    // le plus precis d'abord : "{0} bus, {1} lignes" avant "{0} bus"
    gabarits.sort((a, b) => b.poids - a.poids);
    cache.clear();
  }

  /** L'anglais du texte normalise ``s``, ou null s'il n'est pas connu. */
  function chercher(s, profondeur = 0) {
    if (!LETTRES.test(s)) return null;
    if (cache.has(s)) return cache.get(s);
    let en = exact.has(s) ? exact.get(s) : null;
    // un retour a la ligne du source HTML n'est qu'une espace a l'affichage
    if (en === null && s.includes("\n")) {
      const plat = s.replace(/\n/g, " ");
      if (exact.has(plat)) en = exact.get(plat);
    }
    // les gabarits ne valent que pour des textes courts : un texte de
    // plusieurs milliers de caracteres (code exporte, longue equation) n'est
    // pas un message, et y promener des motifs a trous couterait cher
    if (en === null && s.length <= LONGUEUR_MAX) {
      for (const g of gabarits) {
        if (!s.includes(g.ancre)) continue;
        let m = null;
        try {
          m = g.re.exec(s);
        } catch {
          m = null; // pile du moteur d'expressions depassee : pas de traduction
        }
        if (!m) continue;
        const trous = {};
        g.ordre.forEach((k, i) => {
          const brut = m[i + 1];
          const t = profondeur < 3 ? chercher(normaliser(brut), profondeur + 1) : null;
          trous[k] = t === null ? brut : brut.replace(normaliser(brut), t);
        });
        en = g.en.replace(/\{(\d+)\}/g, (_x, k) => (k in trous ? trous[k] : ""));
        break;
      }
    }
    if (en === null && profondeur < 3 && s.length <= 4 * LONGUEUR_MAX) {
      for (const sep of SEPARATEURS) {
        if (!s.includes(sep)) continue;
        const parts = s.split(sep);
        let connu = false;
        const tr = parts.map((p) => {
          const t = chercher(normaliser(p), profondeur + 1);
          if (t !== null) connu = true;
          return t === null ? p : p.replace(normaliser(p), t);
        });
        if (connu) {
          // " : " francais -> ": " anglais ; les autres separateurs restent
          en = tr.join(sep === " : " ? ": " : sep);
          break;
        }
      }
    }
    if (cache.size > 20000) cache.clear();
    cache.set(s, en);
    if (en === null && profondeur === 0 && inconnus.size < 5000) inconnus.add(s);
    return en;
  }

  /** Le texte dans la langue courante (les espaces de bord sont gardes). */
  function t(texte) {
    if (langue !== "en" || typeof texte !== "string") return texte;
    const cle = normaliser(texte);
    if (!cle) return texte;
    const en = chercher(cle);
    if (en === null) return texte;
    const debut = texte.slice(0, texte.length - texte.trimStart().length);
    const fin = texte.slice(texte.trimEnd().length);
    return debut + en + fin;
  }

  // ------------------------------------------------------------- la page
  const hors = (n) => {
    for (let p = n.nodeType === 1 ? n : n.parentNode; p && p.nodeType === 1; p = p.parentNode) {
      if (HORS.has(p.nodeName.toUpperCase()) || (p.dataset && p.dataset.i18n === "non")) return true;
    }
    return false;
  };

  function traduireTexte(n) {
    if (n.__en !== undefined && n.data === n.__en) return; // notre propre ecriture
    if (hors(n)) return;
    const fr = n.data;
    const en = t(fr);
    n.__fr = fr;
    n.__en = en;
    if (en !== fr) n.data = en;
  }

  function traduireAttributs(e) {
    for (const a of ATTRIBUTS) {
      if (!e.hasAttribute(a)) continue;
      const v = e.getAttribute(a);
      const garde = e.__i18n && e.__i18n[a];
      if (garde && garde.en === v) continue; // notre propre ecriture
      const en = t(v);
      (e.__i18n = e.__i18n || {})[a] = { fr: v, en };
      if (en !== v) e.setAttribute(a, en);
    }
  }

  function traduireSousArbre(racine) {
    if (racine.nodeType === 3) return traduireTexte(racine);
    if (racine.nodeType !== 1 || hors(racine)) return;
    traduireAttributs(racine);
    const w = document.createTreeWalker(racine, NodeFilter.SHOW_TEXT | NodeFilter.SHOW_ELEMENT, {
      acceptNode: (n) =>
        n.nodeType === 1 && (HORS.has(n.nodeName.toUpperCase()) || (n.dataset && n.dataset.i18n === "non"))
          ? NodeFilter.FILTER_REJECT
          : NodeFilter.FILTER_ACCEPT,
    });
    for (let n = w.nextNode(); n; n = w.nextNode()) {
      if (n.nodeType === 3) traduireTexte(n);
      else traduireAttributs(n);
    }
  }

  function rendreFrancais(racine) {
    const w = document.createTreeWalker(racine, NodeFilter.SHOW_TEXT | NodeFilter.SHOW_ELEMENT);
    for (let n = w.nextNode(); n; n = w.nextNode()) {
      if (n.nodeType === 3) {
        if (n.__fr !== undefined && n.data === n.__en) n.data = n.__fr;
        delete n.__fr;
        delete n.__en;
      } else if (n.__i18n) {
        for (const [a, g] of Object.entries(n.__i18n)) if (n.getAttribute(a) === g.en) n.setAttribute(a, g.fr);
        delete n.__i18n;
      }
    }
  }

  function observer() {
    if (observateur) return;
    observateur = new MutationObserver((mutations) => {
      for (const m of mutations) {
        if (m.type === "characterData") traduireTexte(m.target);
        else if (m.type === "attributes") traduireAttributs(m.target);
        else for (const n of m.addedNodes) traduireSousArbre(n);
      }
    });
    observateur.observe(document.documentElement, {
      childList: true, subtree: true, characterData: true, attributes: true, attributeFilter: ATTRIBUTS,
    });
  }

  function appliquer() {
    document.documentElement.lang = langue;
    for (const b of document.querySelectorAll("[data-langue]"))
      b.classList.toggle("actif", b.dataset.langue === langue);
    if (langue === "en") {
      traduireSousArbre(document.documentElement);
      observer();
    } else {
      if (observateur) observateur.disconnect();
      observateur = null;
      rendreFrancais(document.documentElement);
    }
  }

  /** Change la langue de l'interface et la retient pour la prochaine fois. */
  function choisir(l) {
    if (!LANGUES.includes(l) || l === langue) return;
    langue = l;
    try {
      localStorage.setItem(CLE, l);
    } catch {
      // stockage refuse : la langue vaudra pour cette page seulement
    }
    appliquer();
    document.dispatchEvent(new CustomEvent("langue", { detail: l }));
  }

  /** La langue de depart : ?lang= dans l'adresse, sinon celle retenue, sinon
   *  celle du navigateur (francais s'il est en francais, anglais sinon). */
  function initiale() {
    const url = new URLSearchParams(location.search).get("lang");
    if (LANGUES.includes(url)) return url;
    try {
      const gardee = localStorage.getItem(CLE);
      if (LANGUES.includes(gardee)) return gardee;
    } catch {
      // pas de stockage : la langue du navigateur
    }
    return (navigator.language || "fr").toLowerCase().startsWith("fr") ? "fr" : "en";
  }

  // les boites du navigateur ne sont pas dans la page : on les traduit a l'appel
  for (const f of ["alert", "confirm", "prompt"]) {
    const natif = window[f].bind(window);
    window[f] = (message, ...reste) => natif(typeof message === "string" ? t(message) : message, ...reste);
  }

  function demarrer() {
    langue = initiale();
    for (const b of document.querySelectorAll("[data-langue]"))
      b.addEventListener("click", () => choisir(b.dataset.langue));
    appliquer();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", demarrer);
  else demarrer();

  return {
    charger,
    choisir,
    t,
    /** La langue courante : "fr" ou "en". */
    get langue() {
      return langue;
    },
    /** Le parametre de langue des dates et des nombres (Intl). */
    get locale() {
      return langue === "en" ? "en-GB" : "fr-FR";
    },
    /** Le texte francais d'origine d'un noeud ou d'un element traduit. */
    original: (n) => (n.nodeType === 3 ? n.__fr ?? n.data : Array.from(n.childNodes, (c) => I18N.original(c)).join("")),
    /** Les textes montres en anglais sans entree au dictionnaire. */
    manquants: () => Array.from(inconnus),
  };
})();
