# VAR，Local Projectionを用いたインパルス・レスポンスの推定手法とその一致性について(卒論)

Updated at: January 24, 2025 4:51 PM
Tags: shocks, ゼミ

# 1. イントロダクション

財政政策や金融政策の効果を検証する際，政策実施前と実施後のマクロ経済変数の値を比較するだけでは，政策の効果を推定することはできない．なぜならば，マクロ経済変数は長期的に動的かつ相互に影響し合っているため，政策実施前と実施後の単純な比較では政策による純粋な効果のみを取り出すことができないからだ．

そこで，マクロ経済学の実証研究においては，変数間の動的な関係を表現するいくつかのモデルが導入される．モデルには，経済学の理論に基づいて何かしらの制約が課される．この制約により，モデルから “外生的な政策による影響” を取り出すことができる．最終的に，取り出された “政策による外生的な影響” はインパルス・レスポンス（Impulse Response）として可視化される．

本研究はマクロ経済学の実証ツールへの理解を深めることが目的である．まず，モデルとしてベクトル自己回帰モデル（Vector Autoregression; VAR），Local Projectionモデル（LP）を紹介し，識別制約としてリカーシブ制約（recursive restrictions）を解説する．続いて，R言語の lpirfs パッケージに付属している interest_rules_var_data データセットを用いて実際にインパルス・レスポンスをプロットする．最後に条件をそろえた場合，VAR と LP で同様のインパルス・レスポンスが描かれることを解説し実際にプロットして見せる．

# 2. 分析手法

本章では，金融政策の波及効果を例に取り，分析手法について詳しい説明を試みる．

## 2.1 データ

データはR言語の lpirfs パッケージに付属している interest_rules_var_data データセットを利用する．このデータセットには，アメリカ合衆国の1955年第1四半期から2003年第1四半期までのGDPギャップ，インフレ率，金利が含まれている．GDPギャップは実質GDPと潜在GDPのギャップ率，インフレ率はGDPの四半期インフレ率，連鎖加重物価指数（年率％），金利はFF金利の四半期平均（年率％）である．

## 2.2 分析手法

政策の効果は “政策による外生的なショック” を識別することで推定可能である．“外生的” （exogenous）とは経済主体の意思決定とは無関係に発生する変数の変化のことを意味する．“ショック” とは，リーマン・ショックやオイルショックなどマクロ経済変数に予期せぬ変化を与えるものを指す．中央銀行が政策金利を誘導する金融政策も “金融政策ショック”（monetary policy shocks） と呼ばれ，金融政策の効果は “金融政策ショックの波及効果” や “金融政策ショックの伝播”（propagation of shocks）と表現される．また，“識別”（identification）とは，データから興味のあるパラメータを推定することを意味する．特にパラメータの推定が困難な場合に “識別問題” と呼ばれたり，“〜のため識別できない” などというように表現される．

上記で説明した経済学特有の表現は，これ以降注釈なしに使用する．

### ベクトル自己回帰モデル

VAR はSims（1980）によって提唱されたモデルである．まず，一般的なVAR の式を紹介する．

$$
\begin{equation}
Y_t = \sum_{j = 1}^{P} \Phi_j Y_{t-j} + u_t\:\;\;\;
E(u_t)=0\;\;\;\;
E(u_t u_t^\prime)=\Omega 

\end{equation}
$$

$Y_t$ はt期における変数のベクトル，$P$ は$Y_t$を何期前までの$Y$で説明するかを表す “ラグ” ，$\Phi$は回帰係数行列，

$u_t$ はモデルの誤差項ベクトル， $u_t^\prime$ は $u_t$の転置ベクトルである．また，誤差項の期待値は0，誤差項の分散共分散行列を$\Omega$ とする．さらに，回帰係数行列$\Phi_j$は固有値1未満の定常性の仮定が置かれる．

以降，簡単化のためにラグを1，GDPギャップ$y_t$と金利$r_t$ を変数とするVAR(1)の具体例を用いて説明する．VARはラグの数によってVAR(p)のように表現されることがある．

$$
\begin{equation}
\begin{bmatrix}
  y_t \\ r_t
\end{bmatrix} =
\begin{pmatrix}
  \phi_{11} & \phi_{12} \\
  \phi_{21} & \phi_{22}
\end{pmatrix}
\begin{bmatrix}
  y_{t-1} \\ r_{t-1}
\end{bmatrix} +
\begin{bmatrix}
  u_{yt} \\ u_{rt}
\end{bmatrix}
\end{equation}
$$

次の連立方程式は式2と全く同じ式である．VARが複数の変数が将来の変数の値に相互に影響を与え合うモデルであることが理解できる．VARではt-1期とt期の関係を表す回帰係数行列$\Phi$が，t-2期とt-1期や，t-3期とt-2期でも成り立つと仮定している．

$$
\left\{ \,
    \begin{aligned}
    & y_t = \phi_{11}y_{t-1} + \phi_{12}r_{t-1} + u_{yt} \\
    & r_t = \phi_{21}y_{t-1} + \phi_{22}r_{t-1} + u_{rt} 
    \end{aligned}
\right.

$$

式2は現実のデータから例えば最小二乗法（Ordinary Least Square; OLS）などで推定可能な “誘導形” （reduced-form）の式である．しかしながら，後に分かるが、誘導系の誤差項ベクトル $(u_{yt} \;\;u_{yr})^{T}$ からは “金融政策ショックがGDPギャップにどのような影響を及ぼすか” という問いに答えるための有益な情報は得られない．そこで，実際には観察できない真の経済モデルが次の式で与えられると仮定する．

$$
\begin{equation}
\begin{aligned}
\begin{pmatrix}
  \alpha_{11} & \alpha_{12} \\
  \alpha_{21} & \alpha_{22}
\end{pmatrix}
\begin{bmatrix}
  y_t \\ r_t
\end{bmatrix} &=
\begin{pmatrix}
  \beta_{11} & \beta_{12} \\
  \beta_{21} & \beta_{22}
\end{pmatrix}
\begin{bmatrix}
  y_{t-1} \\ r_{t-1}
\end{bmatrix} +
\begin{bmatrix}
\varepsilon^{Demand}_t \\ \varepsilon^{MonPol}_t
\end{bmatrix} \\
E(\varepsilon_t \varepsilon_t^\prime) &= I
\end{aligned}
\end{equation}
$$

式3は，同時点の変数$y_t, r_t$に相関関係があることを仮定している．このように何かしらの構造を仮定するVARを，誘導形と対比して 構造形VAR（Structural VAR; SVAR）と呼ぶ．この式において$(\varepsilon_t^{Demand} \;\;\varepsilon_t^{MonPol})^\prime$ はベクトルの要素がそれぞれ，ディマンドショックと金融政策ショックを表している． $E(\varepsilon_t \varepsilon_t^\prime) = I$ はショックの分散は1，共分散は0とする仮定を表している．式3を変形し，係数行列を再定義すると式4が導ける．

$$
\begin{bmatrix}
  y_t \\ r_t
\end{bmatrix} =
\begin{pmatrix}
  \alpha_{11} & \alpha_{12} \\
  \alpha_{21} & \alpha_{22}
\end{pmatrix}^{-1}
\begin{pmatrix}
  \beta_{11} & \beta_{12} \\
  \beta_{21} & \beta_{22}
\end{pmatrix}
\begin{bmatrix}
  y_{t-1} \\ r_{t-1}
\end{bmatrix} +
\begin{pmatrix}
  \alpha_{11} & \alpha_{12} \\
  \alpha_{21} & \alpha_{22}
\end{pmatrix}^{-1}
\begin{bmatrix}
\varepsilon^{Demand}_t \\ \varepsilon^{MonPol}_t
\end{bmatrix}
$$

$$
\begin{equation}
\begin{bmatrix}
  y_t \\ r_t
\end{bmatrix} =
\begin{pmatrix}
  \phi_{11} & \phi_{12} \\
  \phi_{21} & \phi_{22}
\end{pmatrix}
\begin{bmatrix}
  y_{t-1} \\ r_{t-1}
\end{bmatrix} +
\begin{pmatrix}
  b_{11} & b_{12} \\
  b_{21} & b_{22}
\end{pmatrix}
\begin{bmatrix}
  \varepsilon_{t}^{Demand} \\      
  \varepsilon_{t}^{MonPol}
\end{bmatrix}
\end{equation}
$$

ここで，式2と式4から実際に観測される$(u_{yt} u_{yr})^{T}$についての誘導形の式を導くことができる．

$$
\begin{equation}
\left\{ \,
\begin{aligned}
  u_{yt} &=     
  b_{11}\varepsilon_t^{Demand} +   
  b_{12}\varepsilon_t^{MonPol} \\
  u_{rt} &= 
  b_{21}\varepsilon_t^{Demand} + 
  b_{22}\varepsilon_t^{MonPol}
\end{aligned}
\right .
\end{equation}
$$

この式5の $u_{yt}$を確認すると，誤差項 $u_{yt}$には金融政策ショック $\varepsilon_t^{MonPol}$だけでなくディマンドショック $\varepsilon_t^{Demand}$ も含まれる．したがって，式2の誤差項ベクトルから直接金融政策ショックの効果を推定できないことがわかる．

一方で，この式5からは，金融政策ショックがGDPギャップに与える影響が $b_{12}$で表すことができるとわかる．なぜなら，金融政策ショック $\varepsilon_t^{MonPol}$ は係数 $b_{12}$倍だけ誤差項 $u_{yt}$を増加（もしくは減少）させ，結果的に誤差項を通して $y_t$ の値を変化させるからである．

式5を行列表現に直すと次のようになる．

$$
\begin{equation}
u_t = B \varepsilon_t
\end{equation}
$$

したがって，ショックがマクロ経済変数に与える影響を推定する問題は，誤差項 $u_t$をいかにしてショック $\varepsilon_t$に分解するのかという問題に置き換えられる．

ここで，式6の左辺の誘導形の誤差項と右辺の構造形の関係を利用すると，式1の分散共分散行列を次のように展開できる．

$$
\begin{equation}
\Omega = E[u_t u_t^\prime] = E[B\varepsilon_t\varepsilon_t^\prime B^\prime] = B \Sigma_\varepsilon B^\prime = BIB^\prime = BB^\prime
\end{equation}
$$

式7より，データから直接的に推定できる$(u_{yt}\; u_{rt})$を使って，直接は観測できないはずの構造ショック$(\varepsilon_{yt}\; \varepsilon_{rt})$を識別する問題は， $\Omega = B B^\prime$ を満たす 行列$B$ を探すことに集約された．

しかしながら，これを満たす行列Bは無数に考えられる．誤差項の分散共分散行列を書き下すと，

$$
\begin{equation}
\begin{pmatrix}
  \sigma_y^2 & \sigma_{yr}^2 \\
  \sigma_{yr}^2 & \sigma_r^2
\end{pmatrix} = 
\begin{pmatrix}
  b_{11} & b_{12} \\
  b_{21} & b_{22}
\end{pmatrix}
\begin{pmatrix}
  b_{11} & b_{21} \\
  b_{12} & b_{22}
\end{pmatrix}
\end{equation}
$$

これを連立方程式の形に書き直すと，

$$
\begin{equation}
\left\{ \,
    \begin{aligned}
    & \sigma_y^2 = b_{11}^2 + b_{12}^2 \\
    & \sigma_{yr}^2 = b_{11}b_{21} + b_{12}b_{22} \\
  & \sigma_{yr}^2 = b_{11}b_{21} + b_{12}b_{22} \\
    & \sigma_r^2 = b_{21}^2 + b_{22}^2
    \end{aligned}
\right.
\end{equation}

$$

式9の2つ目と3つ目の式は同じなので，連立方程式は式を3つしか持たない．したがって，行列Bの4つの要素の解を一意に求めることができない．そこで，経済学の理論に基づいて行列 Bに何かしらの仮定を置くことで4つ目の関係式を追加する．この追加される式は “制約” や “識別制約”，もしくは “制約条件” と呼ばれる．

次の項から具体的な識別制約として，リカーシブ制約（recursive restriction）について解説する．

### リカーシブ制約

リカーシブ制約とは，例えば，金融政策ショックはGDPギャップに対して同時的な影響を持たないと仮定することである．つまり，t期の金融政策ショック$\varepsilon_{rt}$は次期のGDPギャップ$y_{t+1}$には影響を与えても，同時期のGDPギャップ$y_t$には影響を与えないということを意味する．金融政策は一般的にその効果が発揮されるまでに時間がかかるとされているため，このような仮定を置く．

一方で，ディマンドショックが金利に同時期的な影響を与えると仮定するのは，今回使用するデータが四半期データなのに対して，中央銀行が金融政策を変更し得るタイミングが4回よりも多く存在するからである．本論文で使用するデータはアメリカのマクロ経済データであるが，アメリカでは金融政策を決定する会合である連邦公開市場委員会（Federal Open Market Committee; FOMC）が年に８回開かれる．また，日本でも同様に金融政策決定会合（Monetary Policy Meetings; MPM）が年に８回開催される．したがって，政府は四半期毎のGDPギャップの値を考慮に入れて同時期のうちに市場金利を操作することが可能である．

この仮定は次の式のように$b_{12} = 0$を意味する．

$$
\begin{equation}
\begin{bmatrix}
  y_t \\ r_t
\end{bmatrix} =
\begin{pmatrix}
  \phi_{11} & \phi_{12} \\
  \phi_{21} & \phi_{22}
\end{pmatrix}
\begin{bmatrix}
  y_{t-1} \\ r_{t-1}
\end{bmatrix} +
\begin{pmatrix}
  b_{11} & 0 \\
  b_{21} & b_{22}
\end{pmatrix}
\begin{bmatrix}
  \varepsilon_{t}^{Demand} \\      
  \varepsilon_{t}^{MonPol}
\end{bmatrix}
\end{equation}
$$

この制約により行列Bは識別可能となる．実際に行列Bの各要素の値を求める．式8にリカーシブ制約を課すと，

$$
\begin{equation}
\begin{pmatrix}
  \sigma_y^2 & \sigma_{yr}^2 \\
  \sigma_{yr}^2 & \sigma_r^2
\end{pmatrix} = 
\begin{pmatrix}
  b_{11} & 0 \\
  b_{21} & b_{22}
\end{pmatrix}
\begin{pmatrix}
  b_{11} & b_{21} \\
  0 & b_{22}
\end{pmatrix}
\end{equation}
$$

式11を連立方程式の形に書き下すと，

$$
\begin{equation}
\left\{ \,
    \begin{aligned}
    & \sigma_y^2 = b_{11}^2 \\
    & \sigma_{yr}^2 = b_{11}b_{21} \\
    & \sigma_r^2 = b_{21}^2 + b_{22}^2
    \end{aligned}
\right.
\end{equation}

$$

式12を解くと，行列Bの要素が全て求まるので，リカーシブ制約を用いてショックの識別が可能となることが確かめられた．

$$
\begin{equation}
\left\{ \,
    \begin{aligned}
    & b_{11} = \sigma_y \\
    & b_{12} = 0 \\
    & b_{21} = \frac{\sigma_{yr}^2}{\sigma_y} \\
    & b_{22} = \sqrt{\sigma_r^2 - \frac{(\sigma_{yr}^2)^2}{\sigma_y^2}}
    \end{aligned}
\right.
\end{equation}
$$

式13を用いると2変数のVARは以下のように表すことができる，

$$
\begin{equation}
\begin{bmatrix}
  y_t \\ r_t
\end{bmatrix} =
\begin{pmatrix}
  \phi_{11} & \phi_{12} \\
  \phi_{21} & \phi_{22}
\end{pmatrix}
\begin{bmatrix}
  y_{t-1} \\ r_{t-1}
\end{bmatrix} +
\begin{pmatrix}
  \sigma_y & 0 \\
  \frac{\sigma_{yr}^2}{\sigma_y} & 
  \sqrt{\sigma_r^2 - \frac{(\sigma_{yr}^2)^2}{\sigma_y^2}}
\end{pmatrix}
\begin{bmatrix}
  \varepsilon_{t}^{Demand} \\      
  \varepsilon_{t}^{MonPol}
\end{bmatrix}
\end{equation}
$$

式11〜14の計算をより一般化した計算手法をコレスキー分解（Cholesky decomposition）と呼ぶ．式7のように，ある行列 $\Omega$ を下三角行列 $B$とその転置行列 $B^\prime$の積に分解する際に用いられる計算手法で， $2\times2$の正方行列だけでなく，より行数と列数の大きい正方行列についても同様に分解可能である．コレスキー分解はR や Python などのプログラミング言語の関数を使用すれば1行で実装することができる．

本項では，２変数のVARを例にリカーシブ制約について説明を試みたが，３変数以上のVARにおいても適切な仮定を増やすことでリカーシブ制約を用いたショックの識別が可能である．

### インパルス・レスポンス

インパルス・レスポンスとは，t期にのみ外生的なショックが発生し，その他の期には一切ショックが発生しなかったと仮定したときに，目的変数 $Y_{t+h}$がt期の外生的なショックから受ける変動のことである．後にわかるが，インパルス・レスポンスはモデルのWold分解表現（もしくはVMA表現）の係数部分そのものである．前項での金融政策ショックを例に挙げれば，t期の外生的なショックは金融政策ショックのみ発生し（ $\varepsilon_t^{Monpol} = 1,\varepsilon_t^{demand}=0$  ），その他の期のショックは一切発生しない（ $\{ \varepsilon_\tau = 0 \}_{t-\infty < t<t,t+1<\tau,t+h}$）と仮定したときに，GDPギャップが $\varepsilon_t^{Monpol}$から受ける変動が分析の対象となるインパルス・レスポンスとなる．インパルス・レスポンスを推定することで，金融政策ショックが発生した時にGDPギャップがどのように変動するのか，その効果はどの程度有意に持続するのか，どれくらいの期間で効果は収束するのか，などを分析できる．

インパルス・レスポンスを求めるには，まず，目的変数 $Y_t$を外生的なショックのベクトル $\{ \varepsilon_\tau \}_{t-\infty < \tau \le t}$ のみで表現する．これは，VARの式を再帰的に代入することで求まる．

$$
\begin{aligned}

Y_t &= \Phi Y_{t-1} + B \varepsilon_t \\

Y_t &= \Phi (\Phi Y_{t-2} + B \varepsilon_{t-1}) + B\varepsilon_t \\

Y_t &= B\varepsilon_t + \Phi B \varepsilon_{t-1} + \Phi^2 Y_{t-2} \\

Y_t &= B\varepsilon_t + \Phi B \varepsilon_{t-1} + \Phi^2 B \varepsilon_{t-2} + \Phi^3 Y_{t-3}　\\

&... \\

Y_t &= B\varepsilon_t + \Phi B \varepsilon_{t-1} + \Phi^2 B \varepsilon_{t-2} + \Phi^3 B \varepsilon_{t-3} + ... + \Phi^\infin Y_{t-\infin} \\

\end{aligned}
$$

$$
\begin{equation}

Y_t = \Phi^\infin Y_{t-\infin} + \sum_{j=1}^{\infin} \Phi^j B \varepsilon_{t-j}

\end{equation}
$$

式15において，この章の初めに行列$\Phi$の固有値は1未満であると仮定しているため $\Phi^\infin$は0に収束する．したがって，$Y_t$はベクトル$\varepsilon_t$のラグ付き変数のみで表現できる．この再帰的に代入することで得られる式をVARモデルのVMA（Vector Moving Average; VMA）表現と呼ぶ．また，より一般的に，ある変数 $Y_t$をそれ以前の期のホワイトノイズのみで表現した式をWold分解表現と呼ぶ．したがって，式16は，VARモデルのVMA表現でありWold分解表現でもある．

$$
\begin{equation}

\begin{aligned}
Y_t &= B\varepsilon_t + \Phi B \varepsilon_{t-1} + \Phi^2 B \varepsilon_{t-2} + \Phi^3 B \varepsilon_{t-3} + ... \\

Y_t &= 
\begin{pmatrix}
b_{11} & b_{12} \\ b_{21} & b_{22}
\end{pmatrix}
\begin{bmatrix}
  \varepsilon_{t}^{Demand} \\      
  \varepsilon_{t}^{MonPol}
\end{bmatrix} \\&+

\begin{pmatrix}
\phi^{I}_{11} & \phi^{I}_{12} \\ 
\phi^{I}_{21} & \phi^{I}_{22}
\end{pmatrix}
\begin{pmatrix}
b_{11} & b_{12} \\ b_{21} & b_{22}
\end{pmatrix} 
\begin{bmatrix}
  \varepsilon_{t-1}^{Demand} \\      
  \varepsilon_{t-1}^{MonPol}
\end{bmatrix}
 \\&+

\begin{pmatrix}
\phi^{I\hspace{-1.2pt}I}_{11} & \phi^{I\hspace{-1.2pt}I}_{12} \\
\phi^{I\hspace{-1.2pt}I}_{21} &
\phi^{I\hspace{-1.2pt}I}_{22}
\end{pmatrix}
\begin{pmatrix}
b_{11} & b_{12} \\ b_{21} & b_{22}
\end{pmatrix} 
\begin{bmatrix}
  \varepsilon_{t-2}^{Demand} \\      
  \varepsilon_{t-2}^{MonPol}
\end{bmatrix}\\&+

\begin{pmatrix}
\phi^{I\hspace{-1.2pt}I\hspace{-1.2pt}I}_{11} & \phi^{I\hspace{-1.2pt}I\hspace{-1.2pt}I}_{12} \\ 
\phi^{I\hspace{-1.2pt}I\hspace{-1.2pt}I}_{21} & \phi^{I\hspace{-1.2pt}I\hspace{-1.2pt}I}_{22}
\end{pmatrix}
\begin{pmatrix}
b_{11} & b_{12} \\ b_{21} & b_{22}
\end{pmatrix}
\begin{bmatrix}
  \varepsilon_{t-3}^{Demand} \\      
  \varepsilon_{t-3}^{MonPol}
\end{bmatrix}

 \\&+ ...

\end{aligned}
\end{equation}
$$

インパルス・レスポンスの値はショックの係数行列で表現されている．式16から金融政策ショックのGDPギャップに対するインパルス・レスポンスは，ホライゾン $h = 0$ のとき$b_{12}$，$h = 1$のとき$\phi^I_{12}b_{22}$，$h = 2$のとき$\phi^{I\hspace{-1.2pt}I}_{12}b_{22}$となり以降も任意のホライゾンまで求めることができる．ホライゾンとはショックが起きた時点から何期先の時点かを表す数のことである．

![スクリーンショット 2025-01-24 15.43.48.png](VAR%EF%BC%8CLocal%20Projection%E3%82%92%E7%94%A8%E3%81%84%E3%81%9F%E3%82%A4%E3%83%B3%E3%83%8F%E3%82%9A%E3%83%AB%E3%82%B9%E3%83%BB%E3%83%AC%E3%82%B9%E3%83%9B%E3%82%9A%E3%83%B3%E3%82%B9%E3%81%AE%E6%8E%A8%E5%AE%9A%E6%89%8B%E6%B3%95%E3%81%A8%E3%81%9D%E3%81%AE%E4%B8%80%E8%87%B4%E6%80%A7%E3%81%AB%E3%81%A4%208b5ba1d7e28a4c0ca98d3d6b4f15a13c/%25E3%2582%25B9%25E3%2582%25AF%25E3%2583%25AA%25E3%2583%25BC%25E3%2583%25B3%25E3%2582%25B7%25E3%2583%25A7%25E3%2583%2583%25E3%2583%2588_2025-01-24_15.43.48.png)

したがって，例えばこのモデルのインパルス・レスポンス関数を $IRF(h)$とするならば，次のように表現できる．

$$
IRF(h) = B\Phi^h
$$

上記のインパルス・レスポンス関数は，金融政策ショックのGDPギャップへのインパルス・レスポンスの他に３通りのインパルス・レスポンスの情報を $2\times2$の行列の成分として持っている．

この項ではVAR(1)のモデルを例に取りインパルス・レスポンスの説明を試みたが，モデルからVMA表現を導くことは任意のラグを持つVARに対して可能である．したがって，どのようなラグ長のVAR(p)についてもインパルス・レスポンスを推定することが可能である．VAR(p)のVMA表現を導く方法については付録を参照されたい．

### Local Projectionモデル

LPはJordà（2005）によって提唱されたモデルである．LPを用いることでVARよりも弱い仮定のもとで直接的にインパルス・レスポンスを求めることができる．

$$
\begin{equation}
Y_{t+h} = \Gamma_h Y_t + \sum_{j = 1}^P A_j Y_{t-j} + u_h
\end{equation}
$$

Yは内生変数のベクトル，$\Gamma_h$は回帰係数行列，$\sum_{j = 1}^p A_j Y_{t-j}$はt期以前の内生変数が$Y_{t+h}$に与える影響を取り除くためのコントロール群である．LPでは，ショックの時点 t期から直接，任意のホライゾン t + h期への影響を推定する．つまり，10期先への影響を知りたければ $Y_{t+10}$を$Y_t$に回帰し，同様に，15期先への影響を知りたければ$Y_{t+15}$を$Y_t$に回帰すれば良い．

例えば$h = 20$ までのショックの波及効果を推定したい場合は，次のように$h = 1$から$h = 20$まで，逐次回帰を行う．

$$
\begin{equation}
\left\{ \,
\begin{aligned}

Y_{t+1} &= \Gamma_1 Y_t + \sum_{j=1}^{P} A_j Y_{t-j} + u_1\\

Y_{t+2} &= \Gamma_2 Y_t + \sum_{j=1}^{P} A_j Y_{t-j} + u_2\\

...\\

Y_{t+20} &= \Gamma_{20} Y_t + \sum_{j=1}^{P} A_j Y_{t-j} + u_{20}\\

\end{aligned}
\right .
\end{equation}
$$

式21の回帰係数行列 $\Gamma_1$から$\Gamma_{20}$はそのままインパルス・レスポンスを表している．ここからは，なぜ$\Gamma_h$がそのままインパルス・レスポンスとなるのかを説明する．

$$
\begin{equation}
\begin{aligned}

Y_{t+10}& = \Gamma_{10} Y_t + u_{10} \\

&= u_{10} + \Gamma_{10} u_{0} + \Gamma_{10}^2Y_{t-10} \\

&...\\

&= \Gamma_{10}^\infin Y_{t-\infin} + u_{10} + \Gamma_{10}u_0 + \Gamma_{10}^2 u_{-10} + ...\\

\end{aligned}
\end{equation}
$$

式22は，コントロール項を省略して$Y_{t + 10}$を$Y_t$に回帰した式について，再帰的に代入を行なった結果である．VARと同様に$\Gamma_h$の固定値が1未満という仮定を置き，誤差項ベクトルをリカーシブ制約のもとで$u_h = B_h \varepsilon_h$に分解可能とすると，式22は次のようになる．

$$
\begin{equation}
Y_{t+10} = IB_{10}\varepsilon_{10} + \Gamma_{10}B_{10}\varepsilon_{0} + \Gamma_{10}^2 B_{10} \varepsilon_{-10} + ...
\end{equation}
$$

ショックはt期に1度だけ発生すると考えるので，式23から$h = 10$の時のインパルス・レスポンスは$\Gamma_{10} B_{10}$と求まることがわかった．式23と同様の式は任意のホライゾン全てで成り立つため，LPにおいて回帰係数行列$\Gamma_h$がそのままインパルス・レスポンスとなることがわかる．

### VARとLPは同一のインパルス・レスポンスを推定している

Plagborg-Møller and Wolf（2021）は，比較的少ない仮定の下でVARとLPが推定するインパルス・レスポンスが一致することを示し，両モデルが最良の推定量を求める線形射影の一種であることを証明した．特にこの主張を無限ラグのVAR($\infin$)だけでなく，有限ラグのVARやリカーシブ制約以外の識別手法を用いる場合にも拡張して示した．本項では前項などで解説したVAR，LP両モデルをより一般化した上で，VAR($\infin$)とLPが同一のインパルス・レスポンスを推定することを確かめる．有限ラグにおける証明は付録を参照されたい．

まず，次のようなデータが観測できるとする．

$$
w_t = (r_t^\prime, x_t, y_t, q_t^\prime)^\prime
$$

$r_t,q_t$はそれぞれ $n_r \times 1$ ， $n_q \times 1$ のコントロール変数のベクトルである． $x_t,y_t$はスカラーである． したがって，$w_t$は $(n_r + 2 + n_q)\times1$のベクトルである．ここでは， $x_t$に生じたインパルスに対する $y_t$の動学的なレスポンスに興味があるとする．続いて，この観測されたデータ$w_t$について次の２つの仮定を課す．

> **仮定１**
データ{ $w_t$ } は共分散定常かつ純粋に非決定的（確率的）で，スペクトル密度行列はどのような周期においても正則行列であり，Wold分解係数は絶対的に総和可能である．
> 

仮定１は，VARとLPの両モデルでインパルス・レスポンスを推定するための最低限の仮定であり，両モデルのインパルス・レスポンスが一致するための追加的な仮定ではない．共分散定常の仮定により，データ $\{w_t\}$は例えば式(15)のようなWold分解した形で表現できる．また，Wold分解係数が総和可能であるとは，係数が無限に発散することなくインパルス・レスポンスは必ず推定できるということである．仮定１についての詳しい説明は付録を参照されたい．

> **仮定２**
データ{$w_t$}は各成分が確率的に正規分布に従う時系列ベクトルである．
> 

仮定2は，条件付き期待値関数（Conditional Expectation Function; CEF）がパラメータについて線形な関数としてかけるようにするための仮定である．一般にCEFの関数型は不明なため，線形射影はあくまでもCEFの近似である．ただし，データが正規分布に従うと仮定するとCEFが線形射影と一致することがわかっている．仮定２についての詳しい説明は付録を参照されたい．

ここで，$w_t$を用いてLPを一般化すると次のようになる．

$$
\begin{equation}
y_{t+h} = \mu_h + \beta_h x_t + \gamma_h^\prime r_t + \sum_{\ell = 1}^\infin \delta_{h,\ell}^\prime w_{t-\ell} + \xi_{h,t}
\end{equation}
$$

式24において，$\xi_{h,t}$は誤差項，$\mu_h, \beta_h, \gamma_{h,1}, \delta_{h,2}, ...$は回帰係数である．$x_t$ は必ずしもショックを表す変数である必要はないことに注意をされたい．また，式24を式20のような形にするには， $r_t$を $x_t$以外のすべての内生変数のベクトルとし $q_t = \mathbf{0}$すれば良い．逆に式20を $\{w_t\}$で表現すると次のように表せる． $i$はベクトルの何番目の要素かを表す添字， $\check\beta$は $(n_r + 2)\times(n_r+2)$の係数行列で $\check{\beta}_{i,\bullet,h}$は 行列 $\beta$から i行目の要素を全て指定している．

$$
w_{i,t+h} = \check{\mu}_{i,h} + \check{\beta}_{i,\bullet, h}w_{t} + \sum_{\ell = 1}^\infin \check{\delta}^\prime_{h,\ell}w_{t-\ell} + \check{\xi}_{h,t}
$$

前項で確認したようにLPにおいては回帰係数がそのままインパルス・レスポンスとなるため，LPインパルス・レスポンスは次のように定義することができる．

> **定義１
$x_t$** に対する ****$y_t$のLPインパルス・レスポンス関数は 式24の $\{\beta_h\}_{h \geq 0}$ で与えられる
> 

したがって，ホライゾンhにおけるLPインパルス・レスポンスを次のように表現することができる．

$$
\begin{equation}
\beta_h = E(y_{t+h} | x_t = 1, r_t, \{w_\tau \}_{\tau < t}) -

E(y_{t+h}|x_t = 0, r_t, \{ w_\tau \}_{\tau<t})
\end{equation}
$$

この式において，$r_t$の同時期の値はコントロール変数に含まれるが，$q_t$の同時期の値は含まれない．コントロール変数を $r_t, q_t$ に分類した理由はこのためである．

続いて，VAR(∞)の一般形を$w_t$を用いて定義すると次のようになる．

$$
\begin{equation}
w_t = c + \sum_{t-1}^\infin A_t w_{t-\ell} + u_t
\end{equation}
$$

ここで，$u_t \equiv w_t - E(w_t | \{ w_\tau \}_{-\infin < \tau < t})$ は誤差項， $A_1, A_2, ...$は射影係数．

誤差項の分散共分散行列をコレスキー分解すると$\Sigma_u = BB^\prime$と表せ，Bは下三角行列となる．式26に対応する再帰的なSVARの表現は次のようになる．

$$
A(L) w_t = c + B \eta_t
$$

$A(L)\equiv I - \sum_{\ell = 1}^\infin A_\ell L^\ell$ で，$\eta_t \equiv B^{-1}u_t$である．$r_t$がVARの最初に順序づけられ、$q_t$は最後に順序づけられることに注意されたい．ラグ多項式を定義する．$\sum_{\ell=0}^\infin C_\ell L^\ell = C(L) = A(L)^{-1}$．$x_t$と$y_t$は $w_t$の$(n_r+1)$番目と$(n_r+2)$番目の要素である（ $w_t$が $(n_r + 2 + n_q)\times 1$のベクトルだったことを思い出してほしい）．ここで，VAR(∞)のインパルス・レスポンス関数を定義すると

> **定義２**
$x_t$のイノベーションに関する$y_t$のVARインパルス・レスポンス関数は$\{ \theta_h \}_{h \ge0}$で与えられる．
$\theta_h \equiv C_{n_r+2, \bullet,h}B_{\bullet, n_r+1},$
> 

![スクリーンショット 2025-01-24 15.01.54.png](VAR%EF%BC%8CLocal%20Projection%E3%82%92%E7%94%A8%E3%81%84%E3%81%9F%E3%82%A4%E3%83%B3%E3%83%8F%E3%82%9A%E3%83%AB%E3%82%B9%E3%83%BB%E3%83%AC%E3%82%B9%E3%83%9B%E3%82%9A%E3%83%B3%E3%82%B9%E3%81%AE%E6%8E%A8%E5%AE%9A%E6%89%8B%E6%B3%95%E3%81%A8%E3%81%9D%E3%81%AE%E4%B8%80%E8%87%B4%E6%80%A7%E3%81%AB%E3%81%A4%208b5ba1d7e28a4c0ca98d3d6b4f15a13c/%25E3%2582%25B9%25E3%2582%25AF%25E3%2583%25AA%25E3%2583%25BC%25E3%2583%25B3%25E3%2582%25B7%25E3%2583%25A7%25E3%2583%2583%25E3%2583%2588_2025-01-24_15.01.54.png)

$C_{i,\bullet,h}$は$C_h$のi行目，$B_{\bullet,j}$はBのj列目を表す．つまり，定義２の右辺はインパルスレスポンスを表す行列から特定のインパルス・レスポンス， $\eta_{n_r+1,t}$が $y_{t+h}$に与える影響 $\theta_{h}$，を示す要素を抜き出している．

これまでの２つの定義において，LPでもVARでも無限ラグ長のモデルを考えている．有限ラグについての命題と証明は付録を参照されたい．LPとVARのアプローチはこれまでの文献において区別されてきたが，実際のところ同様の推定量を得る（Plagborg-Møller and Wolf 2021）．

> **命題１**
仮定1と2の下では，仮定1と2の下では，LPとVARのインパルス・レスポンス関数は比例する．

$\theta_h = \sqrt{E(\tilde{x}_t^2)}\times\beta_h \;for\; all\; h = 0,1,2,...,where \;\tilde{x_t}\equiv x_t - E(x_t|r_t, \{ w_r \}_{-\infin<\tau<t})$
> 

この命題はつまり，LPインパルス・レスポンス関数は適切に順序づけされた再帰的なVARインパルス・レスポンスとして得られる．同様に，どんな再帰的なVARインパルス・レスポンス関数もまた，適切なコントロール変数のあるLPを通して得られることを意味する．より簡単に言い換えれば，VARとLPのインパルス・レスポンスが一致するということである．

これから命題１の証明を試みる．

まず，Frisch-Waughの定理 により、LPのインパルス・レスポンス$\beta_h$について次が成り立つ．Frisch-Waughの定理の詳しい内容については付録を参照されたい．

$$
\begin{equation}
\beta_h = \dfrac{Cov(y_{t+h}, \tilde{x_t})}{E(\tilde{x_t}^2)}
\end{equation}
$$

続いて，式26のVARは以下のようにWold分解することができる．

$$
w_t = \chi + C(L)u_t = \chi +  \sum_{t-0}^\infin C_t B \eta_t, \;\;\chi \equiv C(1)c.
$$

VARのWold分解表現から $y_{t+h}$のWold分解表現を求めると，

$$
\begin{aligned}
y_{t+h} &= w_{n_r+2,t+h} \\ &= \chi_{n_r+2} + \cdots + C_{n_r+2,\bullet,h} B\eta_{t} + \cdots \\
&= \chi_{n_r+2} + \cdots +  (C_{n_r+2,\bullet,h}B_{\bullet,1}\eta_{1,t} + \cdots + \underline{C_{n_r+2,\bullet,h}B_{\bullet,n_r+1}\eta_{n_r+1,t}} + \cdots + C_{n_r+2,\bullet,h}B_{\bullet,n_r+n_q+2}\eta_{n_r+n_q+2,t}) + \cdots
\end{aligned}
$$

この式を $Cov(y_{t+h}, \eta_{n_r+1,t})$に代入すると， $\eta_{n_r+1}$は分散が１でその他のショックとの共分散が０なので $\eta_{n_r+1,t}$の係数のみが残り，次の式が成り立つ．

$$
\begin{equation}
\theta_h = C_{n_r+2,\bullet,h}B_{\bullet, n_r +1} = Cov(y_{t+h}, \eta_{n_r+1,t})
\end{equation}
$$

続いて，$u_t = B \eta_t$とコレスキー分解の性質から次の式を導出する．

$$
\begin{equation}
\eta_{n_r+1,t} = \dfrac{1}{\sqrt{E(\tilde{u}_{n_r+1,t}^2)}} \times \tilde{u}_{n_r+1,t},
\end{equation}
$$

まず，誤差項ベクトル $u_t$のコレスキー分解を用いて $x$についての誤差項 $u_{n_r+1,t}$ をショックの式に分解する．

$$
\begin{aligned}u_{n_r+1,t} &= B_{n_r+1,\bullet}\eta_t \\ &= B_{n_r+1,1}\eta_{1,t} + \cdots  + B_{n_r+1,n_r+1}\eta_{n_r+1,t} + B_{n_r+1,n_r+2}\eta_{n_r+2,t} + \cdots + B_{n_r+1, n_r+n_q+2}\eta_{n_r+n_q+2,t} \end{aligned}
$$

ここで， $B$はコレスキー分解の性質から下三角行列となるため， $u_{n_r+1,t}$は コントロール変数ベクトル$r$についてのショックの部分 $B_{n_r+1,1:n_r}\eta_{1:n_r,t}$と $x$についてのショックの部分 $B_{n_r+1,n_r+1}\eta_{n_r+1,t}$に分解される．ただし， $B_{n_r+1,1:n_r}$は行列 $B$の $n_r+1$行目の $1〜(n_r+1)$列目の要素で構成される横ベクトルである．

$$
u_{n_r+1,t} = B_{n_r+1,1:n_r}\eta_{1:n_r,t} + B_{n_r+1,n_r+1}\eta_{n_r+1,t}
$$

また， $B_{n_r+1,1:n_r}\eta_{1:n_r,t}$は， $u_{n_r+1,t}$の $\eta_{1:n_r,t}$のもとでのCEFとして表せる．このCEFは $u_{n_r+1,t}$ の $u_{1:n_r,t}$のもとでのCEFと同じである．したがって，次のように式を変形できる．

$$
\begin{aligned} &u_{n_r+1,t} = E(u_{n_r+1,t}| \eta_{1:n_r,t}) + B_{n_r+1,n_r+1}\eta_{n_r+1,t} \\ \Rightarrow &B_{n_r+1,n_r+1}\eta_{n_r+1,t} = u_{n_r+1,t} - E(u_{n_r+1,t}| \eta_{1:n_r,t}) \\ \Rightarrow &B_{n_r+1,n_r+1}\eta_{n_r+1,t} = u_{n_r+1,t} - E(u_{n_r+1,t}| u_{1:n_r,t}) \\ \Rightarrow &B_{n_r+1,n_r+1}\eta_{n_r+1,t} = \tilde{u}_{n_r+1,t} \end{aligned}
$$

$E(u_{n_r+1,t}| u_{1:n_r,t})$は， $r$の誤差項ベクトル $u_{1:n_r,t}$ を条件とする$x$の誤差項 $u_{n_r+1,t}$のCEFである．これは， $u_{n_r+1,t}$のうち $u_{1:n_r,t}$で説明できる部分を表している．また， $y$と $q$についての誤差項に係る係数が０であるため， $E(u_{n_r+1,t} | u_{1:n_r,t})　= E(u_{n_r+1,t} | u_{1:n_r,t},\;u_{n_r+2:n_r+n_q+2,t})$が成り立つ．したがって， $\tilde{u}_{n_r+1,t}$は $x$以外の誤差項で説明できない残差であり，この残差は $B_{n_r+1,n_r+1}\eta_{n_r+1,t}$と等しい．

ショックは定義から， $E(\eta_{n_r+1,t}^2) = 1$より，

$$
\begin{aligned}&B_{n_r+1,n_r+1}\eta_{n_r+1,t} = \tilde{u}_{n_r+1,t} \\ \Rightarrow &B^2_{n_r+1,n_r+1}\eta^2_{n_r+1,t} = \tilde{u}^2_{n_r+1,t} \\ \Rightarrow &B^2_{n_r+1,n_r+1}\underbrace{E(\eta^2_{n_r+1,t})}_{\text{=1}} = E(\tilde{u}^2_{n_r+1,t}) \\ \Rightarrow &B_{n_r+1, n_r+1} = \sqrt{E(\tilde{u}_{n_r+1,t}^2)} \end{aligned}
$$

したがって， $\sqrt{E(\tilde{u}_{n_r+1,t}^2)}\eta_{n_r+1,t} = \tilde{u}_{n_r+1,t}$より，式29が求まった．

最後に，次の式を導出する．

$$
\begin{equation}
\tilde{u}_{n_r+1,t}\equiv u_{n_r+1,t} - E(u_{n_r+1,t} | u_{1:n_r,t}) = \tilde{x_t}.
\end{equation}
$$

まず， $\tilde{x}_t \equiv x_t - E(x_t | r_t, \{ w_\tau\}_{-\infin < \tau < t})$ と $u_{n_r+1,t} \equiv x_t - E(x_t | \{ w_\tau \}_{-\infin < \tau < t})$ より，

$$
\begin{aligned}
u_{n_r+1,t} - \tilde{x}_t &= E(x_t | r_t, \{ w_\tau\}_{-\infin < \tau < t}) - E(x_t | \{ w_\tau \}_{-\infin < \tau < t}) \\ &= E(x_t | r_t, \{ w_\tau\}_{-\infin < \tau < t}) - E(E(x_t | \{ w_\tau \}_{-\infin < \tau < t})|r_t, \{ w_\tau \}_{-\infin < \tau < t}) \\ &= E(x_t- E(x_t|\{w_\tau\}_{-\infin<\tau<t}) | r_t, \{ w_\tau \}_{-\infin < \tau < t}) \\ &= E(u_{n_r+1,t}| r_t, \{w_\tau\}_{–\infin < \tau < t})
\end{aligned}
$$

$r_t$は同時期の変数に対して誤差項と通して影響を与えるので， $E(u_{n_r+1,t}| r_t, \{w_\tau\}_{–\infin < \tau < t}) = E(u_{n_r+1,t}| u_{1:n_r,t})$より，次のように式30が導出できる．

$$
\begin{aligned}
u_{n_r+1,t} - \tilde{x}_t &= E(u_{n_r+1,t}|u_{1:n_r,t}) \\ u_{n_r+1,t} - E(u_{n_r+1,t}|u_{1:n_r,t}) &= \tilde{x}_t \\ \tilde{u}_{n_r+1,t} &= \tilde{x}_t 
\end{aligned}
$$

式28,29,30から，次のようにまとめらる．

$$
\theta_h = \dfrac{Cov(y_{t+h}, \tilde{x}_t)}{\sqrt{E(\tilde{x}_t^2)}}
$$

したがって，命題は証明された．

# 3. 結果

次のグラフは，R言語の lpirfs パッケージに付属している interest_rules_var_data データセットを用いて，アメリカの1995年第一四半期から2003年第一四半期までのGDPギャップ，インフレ率，FFレートの３変数のデータを用いてにVAR(2)とVAR(12)におけるインパルス・レスポンスを推定しプロットしたものである．

![図１：Impulse response to a shock of federal funds rate in VAR(2) and VAR(12)](VAR%EF%BC%8CLocal%20Projection%E3%82%92%E7%94%A8%E3%81%84%E3%81%9F%E3%82%A4%E3%83%B3%E3%83%8F%E3%82%9A%E3%83%AB%E3%82%B9%E3%83%BB%E3%83%AC%E3%82%B9%E3%83%9B%E3%82%9A%E3%83%B3%E3%82%B9%E3%81%AE%E6%8E%A8%E5%AE%9A%E6%89%8B%E6%B3%95%E3%81%A8%E3%81%9D%E3%81%AE%E4%B8%80%E8%87%B4%E6%80%A7%E3%81%AB%E3%81%A4%208b5ba1d7e28a4c0ca98d3d6b4f15a13c/Untitled.png)

図１：Impulse response to a shock of federal funds rate in VAR(2) and VAR(12)

$h = 24$までのVAR(12)のインパルス・レスポンスを確認すると，政策金利の上昇（金融政策ショック）に対するGDPギャップのインパルス・レスポンスは２期目には大きくに減少方向に反応し，その後10期先まで減少したのち16期目のあたりで影響が収束していることがわかる．これはつまり，景気は金利の上昇に対して半年という比較的早く反応し，その影響は2年半後に最も大きくなり，大体4年後には影響が金利上昇の影響がなくなるという事である．それに対して，インフレ率のインパルス・レスポンスは6期目から減少方向に反応し，10期目に最も反応が大きくなり，だんだんと反応は小さくなるが24期目でも多少金利上昇の影響が残っている．これはつまり，物価は1年半後にようやく金利上昇に反応し始め，2年半後に最も影響が大きく現れ，その後だんだんと影響が薄れていると解釈できる．この分析結果からは例えば次のようなインプリケーションが得られる．物価は金利に対して反応が遅れる傾向があるので，中央銀行が物価を安定させるためには物価上昇を観測するよりも早く政策金利を上げるなどの意思決定が必要となる．

![図２：Impulse response to all shock of all endogenous variables in LP](VAR%EF%BC%8CLocal%20Projection%E3%82%92%E7%94%A8%E3%81%84%E3%81%9F%E3%82%A4%E3%83%B3%E3%83%8F%E3%82%9A%E3%83%AB%E3%82%B9%E3%83%BB%E3%83%AC%E3%82%B9%E3%83%9B%E3%82%9A%E3%83%B3%E3%82%B9%E3%81%AE%E6%8E%A8%E5%AE%9A%E6%89%8B%E6%B3%95%E3%81%A8%E3%81%9D%E3%81%AE%E4%B8%80%E8%87%B4%E6%80%A7%E3%81%AB%E3%81%A4%208b5ba1d7e28a4c0ca98d3d6b4f15a13c/Untitled%201.png)

図２：Impulse response to all shock of all endogenous variables in LP

図２はLPによって推定されたインパルス・レスポンスを全てプロットしたものである．灰色の帯は95%信頼区間を表している．

![図３：Impulse response to all shock of federal funds rate in VAR(12) and LP](VAR%EF%BC%8CLocal%20Projection%E3%82%92%E7%94%A8%E3%81%84%E3%81%9F%E3%82%A4%E3%83%B3%E3%83%8F%E3%82%9A%E3%83%AB%E3%82%B9%E3%83%BB%E3%83%AC%E3%82%B9%E3%83%9B%E3%82%9A%E3%83%B3%E3%82%B9%E3%81%AE%E6%8E%A8%E5%AE%9A%E6%89%8B%E6%B3%95%E3%81%A8%E3%81%9D%E3%81%AE%E4%B8%80%E8%87%B4%E6%80%A7%E3%81%AB%E3%81%A4%208b5ba1d7e28a4c0ca98d3d6b4f15a13c/Untitled%202.png)

図３：Impulse response to all shock of federal funds rate in VAR(12) and LP

図３は，ラグを12とした際のVARとLPのインパルス・レスポンスをプロットしたものである．金融政策ショックに対するGDPギャップのインパルス・レスポンスを確認すると15期目まではVARとLPでほとんと一致していることが確認できる．インフレ率は，10期目までは概ね一致している．アメリカの政策金利であるFFレートを確認すると，VARとLPで3期目から比較的推定値がずれている．Plagborg-Møller and Wolf（2021）は，有限ラグのモデルにおいても $h = p$ 期まではインパルス・レスポンスの推定量が一致することを示したが，これはあくまでも母集団において成り立つ性質のため，実際のデータから推定される推定量には推定量にはずれがある事がわかる． 

# 4. 結論

本研究はマクロ経済学の実証ツールへの理解を深めることを目的とした．まず，モデルとしてベクトル自己回帰モデル（Vector Autoregression; VAR），Local Projectionモデル（LP）を紹介し，識別制約としてリカーシブ制約（zero contemporaneous restrictions）を解説した．続いて，R言語の lpirfs パッケージに付属している interest_rules_var_data データセットを用いて実際にインパルス・レスポンスをプロットした．最後に母集団においてVARとLPのインパルス・レスポンスが一致することを解説し，VAR と LP で概ね同様のインパルス・レスポンスが描かれることを解説し実際にプロットして見せた．

# 参考文献

Jordà, Ò.. 2005.. “Estimation and Inference of Impulse Responses by Local Projections,” *American Economic Review*, 95 1, 161-182. 

Sims CA. 1980.. “Macroeconomics and reality,” *Econometrica,* 48 1, 1–48. 

Plagborg-Møller, M., and Wolf, CK. 2021.. “Local Projection and VARs Estimate the Same Impulse Responses” *Econometrica*, 89, 2, 955-980. 

Angrist, J. D and J. -S. Pischke.  (2013) 『「ほとんど無害な」計量経済学ー応用経済学のための実証分析ガイド』大森義明，小原美紀，田中隆一，野口晴子（訳）NTT出版 （原著：“Mostly Harmless Econometrics: An Empiricist’s Cimpanion,” *Princeton University Press，* 2008）

北川源四郎（2020），『Rによる時系列モデリング入門』岩波書店．

# 付録

## VAR(p)のVMA表現について

行列表現を用いてVAR(p)をVAR(1)として表現することで，VAR(1)とほとんど同様に再帰的な代入を用いてモデルのVMA表現を得ることができる（Jordà，2005，式5）．

VAR(p)の式を再掲する．

$$
Y_t = \sum_{j = 1}^{P} \Phi_j Y_{t-j} + u_t\:\;\;\;
E(u_t)=0\;\;\;\;
E(u_t u_t^\prime)=\Omega 
$$

この式を次に新しく定義する行列を用いて表す．

$$
W_t^\prime = \begin{bmatrix}Y_t & Y_{t-1} & \cdots & Y_{t-P+1}\end{bmatrix}, \\ F = \begin{pmatrix} \Phi_1 & \Phi_2 & \cdots & \Phi_{P-1} & \Phi_P \\ I & \mathbf{0} & \cdots & \mathbf{0} & \mathbf{0} \\ \mathbf{0} & I & \cdots & \mathbf{0} & \mathbf{0} \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ \mathbf{0} & \mathbf{0} & \cdots & I & \mathbf{0} \end{pmatrix}, \\ U_t^\prime = \begin{bmatrix} B\varepsilon_t & \mathbf{0} & \cdots & \mathbf{0}\end{bmatrix}
$$

$W_t$は $t$期から $t-P+1$期までの内生変数ベクトルの $P$長のベクトルである． $F$は $W_t$を $W_{t-1}$にうまく回帰するための $P\times P$の係数行列である． $U_t$は誤差項ベクトルとゼロベクトルのみで構成される $P$長のベクトルである．これらを用いるとVAR(p)は次のように定義できる．

$$
W_t = FW_{t-1} + U_t
$$

よって，VAR(p)をVAR(1)に変形することができた．この式のVMA表現は次のようになる．

$$
W_t = U_t + FU_{t-1} + F^2 U_{t-2} + \cdots + F^{\infty}U_{t-\infty}
$$

ここで， $W_t$の一番初めの要素についてのみ取り出すとVAR(p)のWold分解表現が得られた．

$$
Y_t = B\varepsilon_t + F_{1,1}B\varepsilon_{t-1} + F_{1,1}^2B\varepsilon_{t-2} + \cdots + F^\infty_{1,1}B\varepsilon_{t-\infty}
$$

## 仮定１の意味について

VAR($\infty$)とLPのインパルス・レスポンスの一致性を説明する際に用いた仮定１について説明する．

仮定１を再掲する．

> **仮定１**
データ{ $w_t$ } は共分散定常かつ純粋に非決定的（確率的）で，スペクトル密度行列はどのような周期においても正則行列であり，Wold分解係数は絶対的に総和可能である．
> 

仮定１は，VARとLPの両モデルでインパルス・レスポンスを推定するための最低限の仮定であり，両モデルのインパルス・レスポンスが一致するための追加的な仮定ではない．まず，データが共分散定常かつ純粋に非決定的（確率的）であるという仮定よりWold分解の定理が成り立つ．Wold分解の定理とは，先述の仮定を敷いたときにデータを無限個のホワイトノイズの線形結合に分解でき，この分解表現は一意に定まる，という定理である．ホワイトノイズとは，期待値が０で分散が時点を問わず一定かつ自己共分散が全て０である時系列である．インパルス・レスポンスは，例えばVARであればVMA表現における無限個のショックの係数部分に現れるため，Wold分解が可能であることは必須の条件となる．

式15を再掲する．

$$

Y_t = \Phi^\infin Y_{t-\infin} + \sum_{j=1}^{\infin} \Phi^j B \varepsilon_{t-j}
$$

次に，Wold分解係数が絶対的に総和可能であるとは，そのままインパルス・レスポンスが総和可能であることを意味する．例えば，式15のようなWold分解表現は，ショックの時期がt期から遠ざかれば遠ざかるほど係数行列の累乗の数が増加するため，係数行列の固有値が１以下でなければ値が無限に発散してしまい計算ができなくなる．したがって，Wold分解係数が絶対的に総和可能であることは，インパルス・レスポンスを推定するために必要となる仮定であることがわかる．

## Frisch-Waugh の定理

K+1変量の重回帰モデルが次の式で与えられるとする．

$$
Y_i = \beta_0 + \beta_1 X_{1,i} + \beta_2 X_{2,i} + \cdots + \beta_K X_{K,i} + u_i
$$

このとき，k番目の回帰係数 $\beta_k$は次の式から求められる．

$$
\beta_k = \dfrac{Cov(Y_i,\tilde{X}_{k,i})}{Var(\tilde{X}_{k,i})}
$$

ただし， $\tilde{X}_{k,i}$は $X_{k,i}$をその他の全ての説明変数に回帰した時の残差である．これが，Frisch-Waughの定理である．

この定理を証明する．まず， $\tilde{X}_{k,i}$を定式化する． $\hat{X}_{k,i}$はその他のすべての説明変数で回帰した時の $X_{k,i}$の推定値である．

$$
\tilde{X}_{k,i} = X_{k,i} -\hat{X}_{k,i}, \\ \hat{X}_{k,i} = \gamma_0 + \gamma_1 X_{1,i} + \cdots + \gamma_{k-1} X_{k-1,i} + \gamma_k X_{k+1,i} + \cdots + \gamma_{K-1}X_{K,i}
$$

これらを用いて， $Cov(Y_i, \tilde{X}_{k,i})$を分解する．

$$
Cov(Y_i,\tilde{X}_{k,i}) = Cov(\beta_0 + \beta_1 X_{1,i} + \beta_2 X_{2,i} + \cdots + \beta_K X_{K,i} + u_i,\tilde{X}_{k,i})
$$

ここで， $\tilde{X}_{k,i}$は残差であるため， $E(\tilde{X}_{k,i}) = 0$かつ $X_{k,i}$以外の説明変数と無相関である．したがって共分散がゼロになることから次の式が成り立つ．

$$
Cov(Y_{i}, \tilde{X}_{k,i}) = Cov(\beta_0 + \beta_{k}X_{k,i} + u_i, \tilde{X}_{k,i})
$$

また，重回帰モデルの誤差項である $u_i$はすべての説明変数と無相関であり， $\beta_0$は定数であるから，

$$
\begin{aligned} Cov(Y_i, \tilde{X}_{k,i}) = Cov(\beta_k X_{k,i}, \tilde{X}_{k,i}) \\ \end{aligned}
$$

$\tilde{X}_{k,i} = X_{k,i} - E(X_{k,i}| X_{-k})$より（ $X_{-k}$はk番目以外の説明変数），さらに次のように式変形ができ，Frisch-Woughの定理は証明された．

$$
\begin{aligned} Cov(Y_i,\tilde{X}_{k,i}) &= Cov(\beta_k(E(X_{k,i}| X_{-k}) + \tilde{X}_{k,i}), \tilde{X}_{k,i}) \\ &= \beta_k Cov(E(X_{k,i} | X_{-k}),\tilde{X}_{k,i}) + \beta_k Cov(\tilde{X}_{k,i},\tilde{X}_{k,i}) \\ &= 0 + \beta_k Var(\tilde{X}_{k,i}) \\ \Rightarrow \beta_k &= \frac{Cov(Y_{i}, \tilde{X}_{k,i})}{Var(\tilde{X}_{k,i})}\end{aligned}
$$

また，この定理は Regression Anatomyとも呼ばれる（Angrist and Pischke，2013）．

## CEFと線形射影について

## ラグが有限の場合のインパルス・レスポンスの一致性

次に，ラグが有限の場合のLPとVARのインパルス・レスポンスについて考える．

LPインパルス・レスポンス推定量$\beta_h(p)$は，無限和がラグpで切り捨てられることを除いて，式(24)と同様に投影における$x_t$の係数として定義される．ここでも，すべての係数と残差は最小二乗法による線形射影の結果であると解釈する．

また，ラグpで切り捨てられた無限和を除くVARの式(26)を考える．ここで，

$A_\ell(p), \ell = 1,2,...,p$ と $\Sigma_u(p)$は対応する回帰係数と残差分散を表すとする．$A(L;p) \equiv I - \sum_{\ell=1}^p A_\ell (p)$ とコレスキー分解の$\Sigma_u(p)=B(p)B(p)^\prime$を定義する．さらにラグ多項式を定義すると次のようになる．

$$

\sum_{\ell = 0}^{\infin}C(p)L^\ell = C(L; p)\equiv A(L; p)^{-1}
$$

これより，

$$
\theta_h(p) \equiv C_{n_r + 2,\bullet,h}(p)B_{\bullet,n_r+1}(p);
$$

射影残差を次のように定義する

$$
\tilde{x}_t(P) \equiv x_t - E(x_t| r_t, \{ w_\tau \}_{t-p\le\tau<t}) = x_t - \sum_{\ell=0}^{p} Q_\ell(p)^\prime w_{t-\ell}
$$

(最後の $Q_0(p)$ の$n_q + 2$要素はゼロ)

さらに $Cov^p(.,.)$ 演算子はデータが実際に上記で定義したパラメータ(A(L; p) Σu (p))を持つVAR(p)モデルに従った場合に仮想的に得られるVAR内の任意の変数間の共分散を示すとする．

> **命題2**
仮定1,2を課す．非負整数h, p が $h \le p$ を満たすとする．
すると，
 $\theta_h(p) = \sqrt{E(\tilde{x}_t(p)^2)} \times \beta_h(p) + \phi_h(p),$ 
余りは，
$\phi_h = \{ E(\tilde{x}_t(p)^2) \}^{-1/2} \sum_{\ell = p - h +1}^p \{ Cov(y_{t+h}, w_{t-\ell}) - Cov^p(y_{t+h}, w_{t-\ell}) \}Q_\ell(p)$
> 

つまり，$\phi_h(p)$が0の時，VARとLPのインパルス・レスポンスは一致する．これは，$h \le p$ が成り立つホライゾンにおいて，VARとLPのインパルス・レスポンスが一致することを意味する．