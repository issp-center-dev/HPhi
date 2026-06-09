---
date: 2026-06-09
datetime: 2026-06-09 21:01 JST
model: Claude Opus 4.8 (1M context)
status: review
topic: feature/nbody-spin-canonical (canonical Spin NBodyInterAll/NBodyG, on top of PR #229 head e35e6c28)
summary: |
  canonical Spin への NBodyInterAll/NBodyG 対応の徹底レビュー。
  ソース5本(nbody_interall.c / nbody_correlation.c / makeHam.c / mltplySpin.c / readdef.c)と
  追加テスト3本を確認。論理・MPI・OpenMP・FullDiag いずれも整合し、致命的バグは検出されず。
  serial Lanczos==FullDiag、serial==MPI(np=4, 単/二ランクフリップ・3体)、
  NBodyInterAll==標準InterAll 基底エネルギー一致、OpenMP4スレッド再現性、validation拒否を実機確認。
---

# feature/nbody-spin-canonical レビュー

## 最新結論

**致命的バグは見つからなかった。** canonical Spin パスは SpinGC パスの忠実な移植で、
唯一の本質的追加である「ビット列 → canonical index 変換(`GetOffComp`)」と
「固定 Sz を壊す項の弾き(conservation check)」のいずれも論理・実機の両面で正しい。
commit 可。push はプロジェクト規定どおり事前確認。下記「残件・実装時注意」は非ブロッキングの観察事項。

## 検証範囲

- 対象 diff: working tree（未コミット）, base = e35e6c28。
- 対象ソース: `src/nbody_interall.c`, `src/nbody_correlation.c`, `src/makeHam.c`,
  `src/mltplySpin.c`, `src/readdef.c`, `src/include/nbody_{interall,correlation}.h`。
- 参照コード: `CheckMPI.c`(Tpow/Nsite規約), `bitcalc.c:214 GetOffComp`,
  `xsetmem.c:203 setmem_large`(list_1buf/v1buf 確保), `diagonalcalc.c:197`,
  `mltplySpin.c`(dispatch), `mltplyMPISpin.c`(既存 MPI spin 規約)。

## 確認した正当性ポイント

1. **モデル分岐**: `iCalcModel == Spin` の分岐が全呼び出し元で正しい。
   - multiply は `mltplyHalfSpin`(canonical, line 261) 経由のみが Spin、
     `mltplyHalfSpinGC`(line 762) は SpinGC。`MultiplyNBodyInterAllSpinGC` 内で
     `iCalcModel==Spin → multiply_nbody_pair_spin` に正しく振り分け。
   - FullDiag は `makeHam.c` の `case Spin`(491-563, 呼び出し 516) と
     `case SpinGC`(338-489, 呼び出し 416) で別個に接続。
   - diagonal は `diagonalcalc.c:197` → `SetDiagonalNBodyInterAllSpinGC` 内で Spin 分岐。

2. **Sz 保存チェック** (`CheckNBodyInterAllSpinConservation` / `CheckNBodyGSpinConservation`):
   `delta_nup = Σ(f[1]-f[3])`(f[1]=spin_out, f[3]=spin_in)で正味アップ数変化を計算し、
   非ゼロを拒否。canonical 化が同一サイト演算子の端点(out/in)を保存するため、
   canonical factor 上での総和は raw 積の Sz 変化に一致。diagonal 項は常に 0 で通過。
   NBodyG 側は `NBodyG_IsZero` 項を正しくスキップ。readdef での実行順
   (Validate→Normalize→**Conservation**→Classify→Hermite)も適切。

3. **ビット規約の一貫性**: `CheckMPI.c:747-764` で SpinGC/Spin とも
   inter-process 用 `Tpow[Nsite]=1` から ×2 で再定義。よって
   `apply_nbody_interall_bits` の `site < Nsite`(intra, `lo & Tpow[site]`)/
   else(inter, `ro & Tpow[site]`=myrank ビット)分岐は両モデルで同一に成立。
   validation の `f[0] >= Nsite` 拒否は readdef 時点(Nsite=グローバル)で走り、
   inter-process サイトを誤って弾かない(MPI 実機で site 7 が通ることを確認)。

4. **canonical index 変換**: 新パスは `apply_*_bits` でビット列を求めた後
   `GetOffComp(list_2_1,list_2_2,...,&j_out)` で局所 index に変換。
   `GetOffComp` は `list_2_1[ia]*list_2_2[ib]==0` で sector 外パターンに FALSE を返す
   (bitcalc.c:245)ため、万一 popcount 不一致が来ても破壊せず安全にスキップ。

5. **MPI ランクフリップ**: `multiply_nbody_pair_spin` / `calc_nbodyg_term_spin` は
   まず `idim_max` を交換して `idim_max_buf` を得てから `list_1`/`tmp_v1`(or vec) を
   `idim_max_buf+1` 受信。canonical はランク毎に次元が異なるため必須で、正しく実装。
   受信バッファ `list_1buf`/`v1buf` は `setmem_large` で `MaxMPI_li(idim_max)+1`
   確保済みで溢れない。origin は involutive(`(r^mask)^mask=r`)で Sendrecv が対称、
   idim=0 ランクでもデッドロック/範囲外なし。Hermite ペアの mask は term0 由来で
   両 term 共通(inter フリップ集合が同一)で正しい。

6. **OpenMP**: 新規 3 つの `#pragma omp parallel for default(none)` の
   data-sharing 句は完全(`gcc-15 -fopenmp -fsyntax-only` がエラーなし)。
   1 term 内では state→j_out が単射のため `tmp_v0[j_out]+=` に競合なし、
   `dam_pr` は reduction。

## 実機検証(すべて pass)

- serial: `lanczos_spin_nbody_interall.sh` → Lanczos==FullDiag, OutputHam に対角/非対角項。
- MPI np=4: `mpi_nbody_interall_spin.sh`(S+_6 S-_7 二ランクフリップ)/
  `mpi_nbodyg_spin.sh` → serial と一致。
- 追加 adversarial(L=8/12, np=4): S+_6 S-_5(intra+inter 単一フリップ)+3体 Sz保存項 →
  serial==MPI(差 ~1e-15)。
- 独立物理クロスチェック: NBodyInterAll で表した S+_0 S-_2 + h.c. の基底エネルギーが
  標準 InterAll で表した同一演算子と一致(-2.7087795191, L=6)。
- OpenMP: gcc-15 -fopenmp ビルドで Lanczos==FullDiag、threads=1/4/8 で基底
  エネルギー再現(差 ~1e-14、reduction 順序)。
- validation: `nbody_interall_validation.sh` / `nbodyg_validation.sh` が
  Sz 非保存項を `delta2Sz=2` で拒否。

## 残件・実装時注意(非ブロッキング)

- **Boost 経路**: `mltplySpinGCBoost` は NBody 関数を呼ばないため、SpinGC+Boost で
  NBodyInterAll を与えると黙って無視される(canonical Spin に Boost 経路は無いので
  本 PR の対象外。PR #229 由来の既存挙動)。気になるなら別途 entry guard を検討。
- **TimeEvolution**: NBodyInterAll は TE を明示拒否、NBodyG(測定)は TE で許可。
  これは意図的で SpinGC と一貫。
- **テストの mpisize**: 追加 MPI テストは np=4 固定(`exact:4`)。CLAUDE.md の教訓どおり
  baseline(`-LE consistency|batching`)では走らない label 構成なので
  np16 の site 分割破綻リスクは無い。L=8 は np=4 で site 分割成立を確認済み。

## レビュー履歴

### 第1回レビュー（2026-06-09 21:01 JST）
本文書。致命的バグなし。serial/MPI/OpenMP/FullDiag/物理クロスチェック/validation を実機確認。
