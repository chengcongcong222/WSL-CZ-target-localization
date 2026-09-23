# R3-C2.2A GPT SYNC

UTC：2026-09-23T08:46:33.049990+00:00

## 三态

- KRAKEN 数据基线：VALIDATED
- Yang 原方法公式恢复：PARTIAL
- Yang 在 E-STD 是否可用：UNDECIDED

## 判定

**C2_2A_PRIMARY_FORMULA_NOT_RECOVERED**

Yang 2015 Eq.(1)–Eq.(9) 全文未取得（closed access；官方 PDF/落地页 403；OpenAlex/Unpaywall 无仓库全文；CORE/仅元数据）。摘要仅提供方法框架。Erratum 提供 g(k_m,z_r)=b_m φ_m(z_r) exp(i k_m δ) 的修正形式及 δ 定义，并指出影响 Eq.(1)(5)(6)(9)，但不恢复主文全部原式。b_m 分解、深度 ambiguity Eq.(6)/(9)、SA 求和核、source-level 处理均为 NOT_RECOVERED。未重构任何 estimator。

## 关键

- Eq.(1)-(9) 主文全文未取得
- Eq.(5) 仅有 Erratum 修正形式 g=b_m φ_m(z_r) exp(i k_m δ)
- b_m / observable / source-level / δ 搜索策略 NOT_RECOVERED
- 未写任何 Yang estimator
- 完成即停止

## 下一步（未执行）

取得主文全文后再做忠实实现；在此之前不得性能判定、不得关闭 Yang。
