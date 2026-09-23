#!/usr/bin/env python3
"""Second-pass PDF retrieval using Crossref official links + mirrors."""
from __future__ import annotations

import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2A"
UA = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)
REF = "https://pubs.aip.org/asa/jasa/article/138/3/1678/680391/"


def try_get(name: str, url: str, headers: dict | None = None) -> None:
    h = {
        "User-Agent": UA,
        "Accept": "application/pdf,text/html,*/*",
        "Referer": REF,
    }
    if headers:
        h.update(headers)
    try:
        req = urllib.request.Request(url, headers=h)
        with urllib.request.urlopen(req, timeout=50) as r:
            body = r.read()
            ctype = r.headers.get("Content-Type", "")
        dest = OUT / f"_dl2_{name}.bin"
        dest.write_bytes(body)
        print(f"{name}: OK {len(body)} ctype={ctype} head={body[:12]!r}")
    except Exception as e:
        print(f"{name}: ERR {e}")


def main() -> None:
    urls = [
        (
            "crossref_main_pdf",
            "https://pubs.aip.org/asa/jasa/article-pdf/138/3/1678/15317235/1678_1_online.pdf",
        ),
        (
            "crossref_err_pdf",
            "https://pubs.aip.org/asa/jasa/article-pdf/144/6/3075/14737039/3075_1_online.pdf",
        ),
        (
            "crossref_main_pdf_http",
            "http://pubs.aip.org/asa/jasa/article-pdf/138/3/1678/15317235/1678_1_online.pdf",
        ),
        (
            "silverchair_main",
            "https://aipp.silverchair-cdn.com/aipp/content_public/journal/jasa/138/3/10.1121_1.4929748/3/1678_1_online.pdf",
        ),
        (
            "s2_pdf",
            "https://pdfs.semanticscholar.org/1db5/d2e264be09bbc9718e1c5b4ecb4765c9e55c.pdf",
        ),
        (
            "core_search",
            "https://api.core.ac.uk/v3/search/works?q=doi:10.1121/1.4929748",
        ),
        (
            "fatcat",
            "https://api.fatcat.wiki/v0/release/lookup?doi=10.1121/1.4929748&expand=files",
        ),
        (
            "scholar_archive",
            "https://scholar.archive.org/search?q=%22Source+depth+estimation+based+on+synthetic+aperture%22",
        ),
        (
            "rg_search",
            "https://www.researchgate.net/search/publication?q=Source%20depth%20estimation%20based%20on%20synthetic%20aperture%20beamfoming",
        ),
    ]
    for name, url in urls:
        try_get(name, url)


if __name__ == "__main__":
    main()
