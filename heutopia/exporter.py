import os
from abc import ABCMeta, abstractmethod
from collections.abc import Sequence

import gspread  # type: ignore
from google.oauth2.service_account import Credentials  # type: ignore
from gspread_dataframe import set_with_dataframe  # type: ignore

from heutopia.batch_result import BatchResult


class AbstractExporter(metaclass=ABCMeta):
    @abstractmethod
    def summary(self) -> str:
        pass

    @abstractmethod
    def dump(self, result: BatchResult) -> None:
        """実行結果をレンダリングする

        Parameters
        ----------
        result : BatchResult
            出力対象の実行結果
        """
        pass

    def dump_src(self, src: str) -> None:
        """ソースコードを出力する

        Parameters
        ----------
        src : str
            ソースコード
        """
        pass


class StdoutExporter(AbstractExporter):
    def __init__(self, obj_columns: Sequence[str] | None = None) -> None:
        """各指標の平均値を標準出力

        Parameters
        ----------
        obj_columns : Sequence[str] | None, optional
            default: []
        """
        self.obj_columns = obj_columns if obj_columns is not None else []

    def summary(self) -> str:
        return "StdoutExporter"

    def dump(self, result: BatchResult) -> None:
        print("")
        print("################ RUN SUMMARY ################")
        print("")

        print("Length: {}".format(len(result.df)))
        print("Mean:")
        print(result.mean.to_json(indent=4))

        for column in self.obj_columns:
            if column not in result.df.columns:
                print("WARNING: {} is not in the result column.".format(column))
            else:
                lst = sorted(zip(result.df[column], result.df["input_filename"]))
                print("Worst cases by {}: {}".format(column, lst[:10]))

        print("")
        print("################ SUMMARY END ################")
        print("")


class LocalFileExporter(AbstractExporter):
    def __init__(self, output_dirpath: str) -> None:
        """ローカルファイルに出力

        Parameters
        ----------
        output_dirpath : str
            ローカルの出力先ディレクトリ
        """
        self.output_dirpath = output_dirpath
        os.makedirs(self.output_dirpath, exist_ok=True)

    def summary(self) -> str:
        return "LocalFileExporter: {}".format(self.output_dirpath)

    def dump(self, result: BatchResult) -> None:
        result.df_with_stat.to_excel(
            os.path.join(self.output_dirpath, "summary.xlsx"), index=False
        )

    def dump_src(self, src: str) -> None:
        with open(os.path.join(self.output_dirpath, "main.cpp.old"), "w") as f:
            f.write(src)


class GSheetsExporter(AbstractExporter):
    # Google Sheets に出力

    def __init__(self, run_id: str, gcp_sa_path: str, sheet_key: str):
        """Google Sheets Exporter

        Parameters
        ----------
        run_id : str
            実行 ID （そのまま Sheet のタブ名になる）
        gcp_sa_path : str
            ローカルに保存した GCP のサービスアカウントの credential JSON ファイルパス
        sheet_key : str
            編集対象の Google Sheets ファイルの key
        """

        self.sheet_key = sheet_key

        scopes = [
            "https://www.googleapis.com/auth/spreadsheets",
            "https://www.googleapis.com/auth/drive",
        ]
        credentials = Credentials.from_service_account_file(gcp_sa_path, scopes=scopes)
        gc = gspread.authorize(credentials)
        gs = gc.open_by_key(sheet_key)
        self.worksheet = gs.add_worksheet(run_id, rows=1, cols=1)

    def summary(self) -> str:
        return "GSheetsExporter: {}".format(self.sheet_key)

    def dump(self, result: BatchResult) -> None:
        set_with_dataframe(
            worksheet=self.worksheet,
            dataframe=result.df_with_stat,
            include_index=False,
            include_column_header=True,
            resize=True,
        )
