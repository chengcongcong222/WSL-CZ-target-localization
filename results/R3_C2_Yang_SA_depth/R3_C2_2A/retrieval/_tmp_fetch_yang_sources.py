#!/usr/bin/env python3
"""Fetch Yang 2015/2018 bibliographic + PDF candidates. No formula invention."""
from __future__ import annotations

import json
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2A"
OUT.mkdir(parents=True, exist_ok=True)
UA = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)
EMAIL = "chengcongcong222@gmail.com"


def get(url: str, timeout: int = 45) -> bytes:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": UA,
            "Accept": "application/pdf,application/json,text/html,text/plain,*/*",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def invert_abstract(inv: dict | None) -> str:
    if not inv:
        return ""
    pos: dict[int, str] = {}
    for w, ids in inv.items():
        for i in ids:
            pos[i] = w
    return " ".join(pos[i] for i in sorted(pos))


def main() -> None:
    log: list[str] = []

    # OpenAlex abstracts
    for doi, name in [
        ("10.1121/1.4929748", "main"),
        ("10.1121/1.5081712", "erratum"),
    ]:
        try:
            raw = get(f"https://api.openalex.org/works/https://doi.org/{doi}")
            d = json.loads(raw)
            abs_text = invert_abstract(d.get("abstract_inverted_index"))
            (OUT / f"_abstract_{name}_openalex.txt").write_text(
                abs_text or "NO_ABSTRACT", encoding="utf-8"
            )
            meta = {
                "title": d.get("title"),
                "open_access": d.get("open_access"),
                "has_fulltext": d.get("has_fulltext"),
                "has_content": d.get("has_content"),
                "locations": [
                    {
                        "pdf_url": loc.get("pdf_url"),
                        "landing": loc.get("landing_page_url"),
                        "is_oa": loc.get("is_oa"),
                    }
                    for loc in d.get("locations") or []
                ],
            }
            (OUT / f"_openalex_{name}_meta.json").write_text(
                json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
            )
            log.append(f"OpenAlex {name}: abs_len={len(abs_text)} oa={meta['open_access']}")
            print(log[-1])
            print(abs_text[:400])
        except Exception as e:
            log.append(f"OpenAlex {name}: ERR {e}")
            print(log[-1])

    # Unpaywall / Crossref / EuropePMC
    endpoints = [
        (
            "unpaywall_main",
            f"https://api.unpaywall.org/v2/10.1121/1.4929748?email={EMAIL}",
        ),
        (
            "unpaywall_err",
            f"https://api.unpaywall.org/v2/10.1121/1.5081712?email={EMAIL}",
        ),
        ("crossref_main", "https://api.crossref.org/works/10.1121/1.4929748"),
        ("crossref_err", "https://api.crossref.org/works/10.1121/1.5081712"),
        (
            "eupmc_main",
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:26428805&format=json&resultType=core",
        ),
        (
            "eupmc_err",
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:30599696&format=json&resultType=core",
        ),
    ]
    for name, url in endpoints:
        try:
            body = get(url)
            (OUT / f"_{name}.json").write_bytes(body)
            log.append(f"{name}: OK {len(body)}")
        except Exception as e:
            log.append(f"{name}: ERR {e}")
        print(log[-1])

    # PDF candidates (record status only; do not invent content)
    pdfs = [
        (
            "hindawi_liang",
            "https://downloads.hindawi.com/journals/mpe/2018/7824671.pdf",
        ),
        (
            "wiley_liang",
            "https://onlinelibrary.wiley.com/doi/pdfdirect/10.1155/2018/7824671",
        ),
        (
            "aip_main_guess1",
            "https://pubs.aip.org/asa/jasa/article-pdf/doi/10.1121/1.4929748/15431211/1678_1_online.pdf",
        ),
        (
            "aip_err_guess1",
            "https://pubs.aip.org/asa/jasa/article-pdf/doi/10.1121/1.5081712/15431212/3075_1_online.pdf",
        ),
    ]
    for name, url in pdfs:
        try:
            body = get(url)
            dest = OUT / f"_dl_{name}.bin"
            dest.write_bytes(body)
            head = body[:8]
            log.append(f"PDF {name}: OK {len(body)} head={head!r}")
        except Exception as e:
            log.append(f"PDF {name}: ERR {e}")
        print(log[-1])

    (OUT / "_retrieval_log.txt").write_text("\n".join(log) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
