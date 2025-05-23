# Notion のページを Export して、GitHub に Markdown ファイルをアップロードする

## 方法

1. Notion のページで、右上の３点アイコンから Export を選択（Markdown & Csv）
2. ダウンロードされたZipファイルを解凍して、フォルダ名・ファイル名をわかりやすく変更
3. GitHubのリモートリポジトリに、フォルダごとアップロード（画像などのパスの問題があるので）
4. Commit Change

## 懸念点

### Notion の数式ブロック

Notion の数式ブロックで書きやすさのために改行を入れていると、きちんと表示されません。

- 完全に改行なし
  - 実際のMarkdown
    ```
    $$Y_t = AY_{t-1} + \varepsilon_t$$
    ```
  - 表示

$$Y_t = AY_{t-1} + \varepsilon_t$$

- 一部改行あり
  - 実際のMarkdown
    ```
    $$
    Y_t = 
    AY_{t-1} + 
    \varepsilon_t
    $$
    ```
  - 表示

$$
Y_t = 
AY_{t-1} + 
\varepsilon_t
$$

- 改行(複数)あり
  - 実際のMarkdown
    ```
    $$
    
    Y_t = 
    AY_{t-1} + 
    \varepsilon_t

    $$
    ```
  - 表示

$$

Y_t = 
AY_{t-1} + 
\varepsilon_t

$$

### \infin ではなく \infty
