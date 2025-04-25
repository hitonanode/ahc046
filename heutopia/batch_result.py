from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class BatchResult:
    # 一括実行結果の `pandas.DataFrame` 管理クラス

    df: pd.DataFrame  # 一括実行結果

    @property
    def mean(self) -> pd.Series:
        return self.df.mean(numeric_only=True).round(4)

    @property
    def df_with_stat(self) -> pd.DataFrame:
        df_mean = self.mean
        df_min = self.df.min(numeric_only=True)
        df_max = self.df.max(numeric_only=True)
        df_mean["input_filename"] = "mean"
        df_min["input_filename"] = "min"
        df_max["input_filename"] = "max"

        return pd.concat([pd.DataFrame([df_mean, df_min, df_max]), self.df]).reindex(
            columns=self.df.columns
        )
