# HEutopia: Heuristic contests optimization tools

短期の AtCoder Heuristic Contest で効率的に実験・パラメータ調整を進めるための Python スクリプト等のテンプレートです。

## はじめにやること

### パッケージのインストール

```bash
uv sync
```

### システムパラメータの設定

`config.yaml` に基本的な設定が集約されている。

## ローカルテストケースの全実行

```bash
make run
```

## Optuna によるパラメータチューニングの実行

```bash
make train
```

## C++ コードにオートフォーマッタをかける

```bash
make format
```

## Python スクリプトにオートフォーマッタをかける

```bash
make pyformat
```
