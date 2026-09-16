#!/usr/bin/env python3
"""Write codemeta.json for an R package without an R installation.

Mirrors the output of the R package `codemeta` (codemeta::write_codemeta()),
except that `fileSize` and `runtimePlatform` are omitted: both describe the
machine that generated the file rather than the package, so would change
between runs.

Authors@R and inst/CITATION are R code. Rather than run R, a small evaluator
handles the subset those files use in practice: strings, numbers, variables,
`meta$Field`, c(), paste(), paste0(), sprintf(), person(), as.person(),
personList(), bibentry(), citEntry(), citHeader() and citFooter(). Anything
else stops with an error, so an unsupported construct is never silently
mis-rendered.

Usage: python3 codemeta.py [package_dir]
"""

import gzip
import json
import os
import re
import sys
import time
import urllib.request

BASE_PACKAGES = {
    "R", "base", "compiler", "datasets", "graphics", "grDevices", "grid",
    "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools",
    "utils",
}

CRAN = {
    "@id": "https://cran.r-project.org",
    "@type": "Organization",
    "name": "Comprehensive R Archive Network (CRAN)",
    "url": "https://cran.r-project.org",
}
BIOC = {
    "@id": "https://www.bioconductor.org",
    "@type": "Organization",
    "name": "Bioconductor",
    "url": "https://www.bioconductor.org",
}

# codemeta's extdata/cran-to-spdx.csv, less the entries without an SPDX id.
CRAN_TO_SPDX = {
    "AGPL-3": "AGPL-3.0", "Apache License version 1.1": "Apache-1.1",
    "Apache License version 2.0": "Apache-2.0", "Artistic-1.0": "Artistic-1.0",
    "Artistic-2.0": "Artistic-2.0", "BSD": "BSD-3-Clause",
    "BSD_2_clause": "BSD-2-Clause", "BSD_3_clause": "BSD-3-Clause",
    "BSL version 1.0": "BSL-1.0", "CC BY 3.0": "CC-BY-3.0",
    "CC BY 4.0": "CC-BY-4.0", "CC BY-NC 3.0": "CC-BY-NC-3.0",
    "CC BY-NC 4.0": "CC-BY-NC-4.0", "CC BY-NC-ND 3.0 US": "CC-BY-NC-ND-3.0",
    "CC BY-NC-ND 4.0": "CC-BY-NC-ND-4.0", "CC BY-NC-SA 3.0": "CC-BY-NC-SA-3.0",
    "CC BY-NC-SA 4.0": "CC-BY-NC-SA-4.0", "CC BY-SA 2.0": "CC-BY-SA-2.0",
    "CC BY-SA 3.0": "CC-BY-SA-3.0", "CC BY-SA 3.0 US": "CC-BY-SA-3.0",
    "CC BY-SA 4.0": "CC-BY-SA-4.0", "CC0": "CC0-1.0",
    "CeCILL version 2": "CECILL-2.0", "CPL version 1.0": "CPL-1.0",
    "EPL version 1.0": "EPL-1.0", "EUPL version 1.1": "EUPL-1.1",
    "FreeBSD": "BSD-2-Clause-FreeBSD", "GPL-2": "GPL-2.0", "GPL-3": "GPL-3.0",
    "LGPL-2": "LGPL-2.0", "LGPL-2.1": "LGPL-2.1", "LGPL-3": "LGPL-3.0",
    "Lucent Public License version 1.02": "LPL-1.02",
    "MPL version 1.0": "MPL-1.0", "MPL version 1.1": "MPL-1.1",
    "MPL version 2.0": "MPL-2.0", "Zlib": "Zlib", "MIT": "MIT",
    "GPL": "GPL-2.0", "GPL 2": "GPL-2.0", "GPL 3": "GPL-3.0",
    "LGPL 2": "LGPL-2.0", "LGPL 2.1": "LGPL-2.1", "LGPL 3": "LGPL-3.0",
    "MPL 1.0": "MPL-1.0", "MPL 1.1": "MPL-1.1", "MPL 2.0": "MPL-2.0",
    "AGPL 3": "AGPL-3.0", "Apache License 1.1": "Apache-1.1",
    "Apache License 2.0": "Apache-2.0", "Artistic 1.0": "Artistic-1.0",
    "Artistic 2.0": "Artistic-2.0",
}

ADDITIONAL_TERMS = [
    "affiliation", "applicationCategory", "applicationSubCategory",
    "copyrightYear", "dateCreated", "dateModified", "downloadUrl", "editor",
    "fileSize", "funder", "identifier", "installUrl", "isAccessibleForFree",
    "isPartOf", "keywords", "memoryRequirements", "operatingSystem",
    "permissions", "processorRequirements", "producer", "provider",
    "publisher", "funding", "relatedLink", "releaseNotes", "sameAs",
    "softwareHelp", "sponsor", "storageRequirements", "supportingData",
    "targetProduct", "contIntegration", "buildInstructions",
    "developmentStatus", "embargoDate", "readme", "issueTracker",
    "referencePublication",
]

BIBTYPES = {
    "Article": "ScholarlyArticle", "Book": "Book", "Booklet": "Book",
    "Inbook": "Chapter", "Incollection": "CreativeWork",
    "Inproceedings": "ScholarlyArticle", "Manual": "SoftwareSourceCode",
    "Mastersthesis": "Thesis", "Misc": "CreativeWork", "Phdthesis": "Thesis",
    "Proceedings": "ScholarlyArticle", "Techreport": "ScholarlyArticle",
    "Unpublished": "CreativeWork",
}


class Unsupported(Exception):
    pass


# --------------------------------------------------------------------------
# DESCRIPTION
# --------------------------------------------------------------------------

def read_dcf(path):
    fields, key = {}, None
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line[:1] in (" ", "\t"):
                if key:
                    fields[key] += "\n" + line
            elif ":" in line:
                key, value = line.split(":", 1)
                fields[key] = value.strip()
    return fields


def clean_str(s):
    if s is None:
        return None
    s = re.sub(r"\s+", " ", s).replace("{", "").replace("}", "").strip()
    return s or None


def parse_deps(value):
    out = []
    for item in (value or "").split(","):
        item = item.strip()
        if not item:
            continue
        m = re.match(r"^([^\s(]+)\s*(?:\(\s*([<>=]+)\s*([^)\s]+)\s*\))?$", item)
        if not m:
            raise Unsupported(f"cannot parse dependency {item!r}")
        version = f"{m.group(2)} {m.group(3)}" if m.group(2) else None
        out.append((m.group(1), version))
    return out


# --------------------------------------------------------------------------
# A tiny evaluator for the R used by Authors@R and inst/CITATION
# --------------------------------------------------------------------------

TOKEN = re.compile(r"""
    (?P<ws>\s+|\#[^\n]*)
  | (?P<str>"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')
  | (?P<num>\d+(?:\.\d*)?(?:[eE][+-]?\d+)?L?)
  | (?P<name>`[^`]+`|[A-Za-z.][A-Za-z0-9._]*)
  | (?P<op><-|[(),=$;])
""", re.VERBOSE)

ESCAPES = {"n": "\n", "t": "\t", "\\": "\\", '"': '"', "'": "'"}


def tokenize(src):
    pos, out = 0, []
    while pos < len(src):
        m = TOKEN.match(src, pos)
        if not m:
            raise Unsupported(f"unexpected character {src[pos]!r}")
        pos = m.end()
        kind = m.lastgroup
        text = m.group(kind)
        if kind == "ws" or text == ";":
            continue
        if kind == "str":
            text = re.sub(r"\\(.)", lambda e: ESCAPES.get(e.group(1), e.group(1)),
                          text[1:-1])
        elif kind == "name":
            text = text.strip("`")
        out.append((kind, text))
    return out


class Parser:
    def __init__(self, src):
        self.toks = tokenize(src)
        self.i = 0

    def peek(self, offset=0):
        j = self.i + offset
        return self.toks[j] if j < len(self.toks) else (None, None)

    def take(self, text=None):
        tok = self.peek()
        if text is not None and tok[1] != text:
            raise Unsupported(f"expected {text!r}, found {tok[1]!r}")
        self.i += 1
        return tok

    def program(self):
        out = []
        while self.peek()[0] is not None:
            out.append(self.statement())
        return out

    def statement(self):
        if self.peek()[0] == "name" and self.peek(1)[1] in ("<-", "="):
            name = self.take()[1]
            self.take()
            return ("assign", name, self.expr())
        return self.expr()

    def expr(self):
        kind, text = self.take()
        if kind == "str":
            node = ("lit", [text])
        elif kind == "num":
            n = float(text.rstrip("L"))
            node = ("lit", [str(int(n)) if n.is_integer() else str(n)])
        elif kind == "name":
            node = ("var", text)
        else:
            raise Unsupported(f"unexpected {text!r}")
        while True:
            if self.peek()[1] == "$":
                self.take()
                node = ("field", node, self.take()[1])
            elif self.peek()[1] == "(" and node[0] == "var":
                self.take()
                args = []
                while self.peek()[1] != ")":
                    if self.peek()[1] == ",":  # empty argument
                        self.take()
                        continue
                    if self.peek()[0] == "name" and self.peek(1)[1] == "=":
                        name = self.take()[1]
                        self.take()
                        args.append((name, self.expr()))
                    else:
                        args.append((None, self.expr()))
                    if self.peek()[1] == ",":
                        self.take()
                self.take(")")
                node = ("call", node[1], args)
            else:
                return node


def as_person_string(s, bibtex=False):
    """as.person() on a string: 'Given Family <email> [roles]', people joined
    by ' and '. In BibTeX style (a bibentry author string), 'Family, Given'
    is also accepted."""
    people = []
    for part in re.split(r"\s+and\s+", s):
        email = re.search(r"<([^>]*)>", part)
        roles = re.search(r"\[([^\]]*)\]", part)
        bare = re.sub(r"<[^>]*>|\[[^\]]*\]|\([^)]*\)", "", part).strip()
        if bibtex and "," in bare:
            family, given = (x.strip() for x in bare.split(",", 1))
            name = given.split() + [family]
        else:
            name = bare.split()
        if not name:
            continue
        people.append({
            "given": name[:-1] or None,
            "family": name[-1] if len(name) > 1 else None,
            "email": [email.group(1)] if email else None,
            "role": [r.strip() for r in roles.group(1).split(",")] if roles else None,
            "comment": None,
        })
    return people


class Evaluator:
    def __init__(self, meta):
        self.env = {}
        self.meta = meta

    def run(self, src):
        """Evaluate a script; return the bibentries its statements yield."""
        entries = []
        for stmt in Parser(src).program():
            value = self.eval(stmt)
            if value and isinstance(value[0], dict) and "bibtype" in value[0]:
                entries.extend(value)
        return entries

    def value(self, src):
        """Evaluate a script; return the value of its last statement."""
        value = []
        for stmt in Parser(src).program():
            value = self.eval(stmt)
        return value

    def strings(self, value):
        if not all(isinstance(v, str) for v in value):
            raise Unsupported("expected a character vector")
        return value

    def eval(self, node):
        kind = node[0]
        if kind == "lit":
            return node[1]
        if kind == "assign":
            self.env[node[1]] = self.eval(node[2])
            return self.env[node[1]]
        if kind == "var":
            if node[1] in self.env:
                return self.env[node[1]]
            if node[1] in ("NULL", "NA"):
                return []
            raise Unsupported(f"unknown variable {node[1]!r}")
        if kind == "field":
            if node[1] != ("var", "meta"):
                raise Unsupported("only meta$Field is supported")
            value = self.meta.get(node[2])
            return [value] if value is not None else []
        return self.call(node[1], node[2])

    def call(self, fn, args):
        positional = [self.eval(a) for n, a in args if n is None]
        named = {n: self.eval(a) for n, a in args if n is not None}

        if fn in ("c", "personList"):
            out = []
            for v in positional + list(named.values()):
                out.extend(v)
            if named and fn == "c":
                return [dict(zip(named, (v[0] for v in named.values())))]
            return out
        if fn in ("paste", "paste0"):
            sep = named.get("sep", [" " if fn == "paste" else ""])[0]
            parts = [self.strings(p) for p in positional]
            n = max((len(p) for p in parts), default=0)
            out = [sep.join(p[i % len(p)] for p in parts if p) for i in range(n)]
            if "collapse" in named:
                out = [named["collapse"][0].join(out)]
            return out
        if fn == "sprintf":
            fmt = self.strings(positional[0])[0]
            vals = iter(v[0] for v in positional[1:])
            return [re.sub(r"%[sd]", lambda _: str(next(vals)), fmt)]
        if fn in ("citHeader", "citFooter"):
            return []
        if fn == "person":
            order = ["given", "family", "middle", "email", "role", "comment"]
            p = dict.fromkeys(order)
            for key, value in zip(order, positional):
                p[key] = value
            p.update({k: v for k, v in named.items() if k in p})
            p["given"] = p["given"] or named.get("first")
            p["family"] = p["family"] or named.get("last")
            p["family"] = p["family"][0] if p["family"] else None
            return [p]
        if fn == "as.person":
            out = []
            for v in positional[0]:
                out.extend(as_person_string(v) if isinstance(v, str) else [v])
            return out
        if fn in ("bibentry", "citEntry"):
            entry = {k: v for k, v in named.items() if v}
            if "author" in entry:
                entry["author"] = [q for p in entry["author"] for q in
                                   (as_person_string(p, bibtex=True)
                                    if isinstance(p, str) else [p])]
            entry["bibtype"] = (named.get("bibtype") or named.get("entry"))[0]
            if positional:
                entry["bibtype"] = positional[0][0]
            return [entry]
        raise Unsupported(f"unsupported function {fn}()")


# --------------------------------------------------------------------------
# codemeta terms
# --------------------------------------------------------------------------

def unbox(values):
    if values is None:
        return None
    return values[0] if len(values) == 1 else list(values)


def person_to_schema(p):
    if p.get("family") is None or not p.get("given"):
        name = (p.get("given") or []) + ([p["family"]] if p.get("family") else [])
        out = {"@type": "Organization", "name": unbox(name)}
    else:
        out = {"@type": "Person", "givenName": unbox(p["given"]),
               "familyName": p["family"]}
    if p.get("email"):
        out["email"] = unbox(p["email"])
    orcid = orcid_of(p.get("comment"))
    if orcid:
        out["@id"] = orcid
    return out


def orcid_of(comment):
    if not comment:
        return None
    values = [v for c in comment for v in (c.values() if isinstance(c, dict) else [c])]
    linked = [v for v in values if "orcid" in v]
    if linked:
        return unbox(linked)
    for c in comment:
        if isinstance(c, dict) and "ORCID" in c:
            oid = c["ORCID"]
            return oid if re.match(r"^https?", oid) else "https://orcid.org/" + oid
    return None


def people_with_roles(people, roles):
    out = []
    for role in roles:
        for p in people:
            has = not p.get("role") if role is None else role in (p.get("role") or [])
            if has:
                out.append(person_to_schema(p))
    return out


class Providers:
    # Providers recorded in the previous codemeta.json are reused, so package
    # lists are downloaded only for new dependencies. A recorded absence (not
    # yet on CRAN/Bioconductor) is trusted only while Version is unchanged, so
    # it is re-checked at the next version bump.
    def __init__(self, previous=None, version=None):
        self.cache = {}
        self.known = {}
        try:
            with open(previous, encoding="utf-8") as f:
                old = json.load(f)
        except (OSError, TypeError, ValueError):
            return
        if isinstance(old, dict):
            entries = [old] + [e for key in ("softwareSuggestions", "softwareRequirements")
                               for e in old.get(key) or [] if isinstance(e, dict)]
            for e in entries:
                provider = next((p for p in (CRAN, BIOC) if p == e.get("provider")), None)
                if e.get("identifier") and (provider or old.get("version") == version):
                    self.known[e["identifier"]] = provider

    def packages(self, url):
        if url not in self.cache:
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url + "/src/contrib/PACKAGES.gz",
                                                timeout=60) as r:
                        text = gzip.decompress(r.read()).decode("utf-8", "replace")
                    break
                except OSError:
                    if attempt == 3:
                        raise
                    time.sleep(10 * 2 ** attempt)
            self.cache[url] = set(re.findall(r"^Package:\s*(\S+)", text, re.M))
        return self.cache[url]

    def guess(self, pkg):
        if pkg in BASE_PACKAGES:
            return None
        if pkg in self.known:
            return self.known[pkg]
        if pkg in self.packages("https://cloud.r-project.org"):
            return CRAN
        if pkg in self.packages("https://www.bioconductor.org/packages/release/bioc"):
            return BIOC
        return None


def format_dep(pkg, version, remotes, providers):
    dep = {"@type": "SoftwareApplication", "identifier": pkg, "name": pkg}
    if version:
        dep["version"] = version
    provider = providers.guess(pkg)
    if provider:
        dep["provider"] = provider
    remote = [r for r in remotes if r.endswith("/" + pkg)]
    if remote:
        dep["sameAs"] = "https://github.com/" + remote[0].replace("github::", "")
    elif provider is CRAN:
        dep["sameAs"] = "https://CRAN.R-project.org/package=" + pkg
    elif provider is BIOC:
        dep["sameAs"] = "https://bioconductor.org/packages/release/bioc/html/" + pkg + ".html"
    return dep


def spdx_license(x):
    first = x.split("|")[0].strip()
    first = re.sub(r"(\w+)\s*\([=>]=\s*(\d+\.*\d*)\)", r"\1 \2", first)
    first = re.sub(r"(.*) \+ file (LICENSE|LICENCE)$", r"\1", first)
    spdx = CRAN_TO_SPDX.get(first)
    return "https://spdx.org/licenses/" + spdx if spdx else first


def citation_to_schema(bib):
    get = lambda k: unbox(bib.get(k)) if bib.get(k) else None
    doi = get("doi")
    doi_url = None
    if doi:
        doi_url = doi if doi.startswith("https://doi.org/") else (
            "https://doi.org/" + doi if doi.startswith("10.") else None)
    authors = people_with_roles(bib.get("author") or [], ["aut", None])
    out = {
        "@type": BIBTYPES.get(bib["bibtype"][:1].upper() + bib["bibtype"][1:].lower()),
        "datePublished": get("year"),
        "author": authors or None,
        "name": get("title"),
        "identifier": doi,
        "url": get("url"),
        "description": get("note"),
        "pagination": get("pages"),
        "@id": doi_url,
        "sameAs": doi_url,
    }
    out = {k: v for k, v in out.items() if v is not None}
    if bib.get("journal"):
        volume = {"@type": ["PublicationVolume", "Periodical"],
                  "volumeNumber": get("volume"), "name": get("journal")}
        issue = {"@type": "PublicationIssue", "issueNumber": get("number"),
                 "datePublished": get("year"),
                 "isPartOf": {k: v for k, v in volume.items() if v is not None}}
        out["isPartOf"] = {k: v for k, v in issue.items() if v is not None}
    return out


def codemeta(root, previous=None):
    d = read_dcf(os.path.join(root, "DESCRIPTION"))
    meta = {k: clean_str(v) for k, v in d.items()}
    pkg = d["Package"]
    version = ".".join(str(int(x)) for x in re.split(r"[.-]", d["Version"]))
    providers = Providers(previous, version)

    cm = {
        "@context": "https://doi.org/10.5063/schema/codemeta-2.0",
        "@type": "SoftwareSourceCode",
        "identifier": pkg,
        "description": clean_str(d.get("Description")),
        "name": f"{pkg}: {clean_str(d.get('Title'))}",
    }
    urls = [u for u in re.split(r"[,\s]+", d.get("URL", "")) if re.match(r"^https?://", u)]
    if urls:
        repo = next((u for u in urls if re.search(r"github\.com|gitlab\.com", u)), urls[0]) \
            if len(urls) > 1 else urls[0]
        if len(urls) > 1:
            cm["relatedLink"] = unbox([u for u in dict.fromkeys(urls) if u != repo])
        cm["codeRepository"] = re.sub(r"#.*$", "", repo)
    if d.get("BugReports"):
        cm["issueTracker"] = d["BugReports"]
    if meta.get("License"):
        cm["license"] = spdx_license(meta["License"])
    # As package_version(): components as integers, separated by dots.
    cm["version"] = version
    cm["programmingLanguage"] = {"@type": "ComputerLanguage", "name": "R",
                                 "url": "https://r-project.org"}
    provider = providers.guess(pkg)
    if provider:
        cm["provider"] = provider

    if "Authors@R" not in d:
        raise Unsupported("DESCRIPTION has no Authors@R field; replace the "
                          "legacy Author/Maintainer fields with Authors@R")
    people = Evaluator(meta).value(d["Authors@R"])
    for term, roles in [("author", ["aut", None]),
                        ("contributor", ["ctb", "com", "dtc", "ths", "trl"]),
                        ("copyrightHolder", ["cph"]), ("funder", ["fnd"]),
                        ("maintainer", ["cre"])]:
        found = people_with_roles(people, roles)
        if found:
            cm[term] = found

    remotes = [r for r in re.split(r"[,\s]+", d.get("Remotes", "")) if r]
    suggests = [format_dep(p, v, remotes, providers)
                for p, v in parse_deps(d.get("Suggests"))]
    requires = [format_dep(p, v, remotes, providers)
                for field in d if field in ("Depends", "Imports")
                for p, v in parse_deps(d[field])]
    if suggests:
        cm["softwareSuggestions"] = suggests
    if d.get("SystemRequirements"):
        requires.append(clean_str(d["SystemRequirements"]))
    if requires:
        cm["softwareRequirements"] = requires

    for term in ADDITIONAL_TERMS:
        value = d.get("X-schema.org-" + term)
        if value:
            cm[term] = unbox([re.sub(r"\s+", "", v) for v in value.split(",")])

    for path in ("inst/CITATION", "CITATION"):
        path = os.path.join(root, path)
        if os.path.exists(path):
            with open(path, encoding="utf-8") as f:
                entries = Evaluator(meta).run(f.read())
            if entries:
                cm["citation"] = [citation_to_schema(b) for b in entries]
            break
    return cm


if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "codemeta.json")
    try:
        cm = codemeta(root, previous=out)
    except Unsupported as e:
        sys.exit(f"codemeta.py: {e}.")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(cm, f, indent=2, ensure_ascii=False)
        f.write("\n")
