# PRIMARY_SOURCE_SET
UTC: 2026-09-23T12:44:30.356064+00:00

身份按 **title/DOI** 判定，不按本地文件名。

| tag | local file | SHA256 | pages | identity (title/DOI) | role |
| --- | --- | --- | --- | --- | --- |
| P1_Yang2014 | `2014.pdf` | `bc389afbf9a77c562775bc4ab051b64b031e7ccbb606f633fb6b1f9059c2e01b` | 13 | Data-based matched-mode source localization for a moving source; DOI 10.1121/1.4863270 | Ref.7 environment + KRAKEN condition |
| P2_Erratum2014 | `2018.pdf` | `e04de4c088ef3e94b531eb0edb7255196a79b9a71a04db7b9d85a18e0e90d33f` | 2 | Erratum to Yang 2014; DOI 10.1121/1.4919288; **490 m -> 4990 m** | span correction |
| P3_Yang2015 | `yang2015.pdf` | `b34952ead13fdaabdf9f1e7d87728fc2cb72fd09497240aac20ebdfc47f4b5b4` | 10 | Source depth estimation based on synthetic aperture beamfoming for a moving source; DOI 10.1121/1.4929748 | core method Eq.(1)-(10), App.A1-A4 |
| P4_Erratum2015 | `2018-2.pdf` | `5123bf0278a55c94e181428a3977247ba7bb9cc3b93d5720238fa8e54873a4da` | 2 | Erratum to Yang 2015; DOI 10.1121/1.5081712; offset-range delta | delta mechanism; paper sim delta=0 |

## 历史状态更新

- `C2_1_PARSER_VALIDATED` / `C2_2B_METHOD_SPEC_LOCKED` 冻结
- `PAPER_REPRO_NEEDS_REF7` -> **SUPERSEDED_BY_USER_REF7**
- `PAPER_REPRO_ENV` -> **AVAILABLE_FROM_PRIMARY_SOURCES**
- `CARRY_FORWARD_FROM_PRIOR_ERRATUM_RECOVERY` -> **ERRATUM_PRIMARY_SOURCE_LOCKED**
- `YANG_ROUTE_UNDECIDED` 保持

## P2 勘误原文（关键）

On page 1222, line 22 in Sec. III C, "490 m" should be "4990 m."

## P4 勘误原文（关键）

delta = r1 - r1_data = r2 - r2_data;
g(k_m,z_r)= b_m phi_m(z_r) exp(i k_m delta);
paper simulations assume delta=0; real data must search initial range / offset-range.
