.PHONY: build format pyformat run train show-studies clean

build: main.cpp
	g++-14 -Wfatal-errors -Wall -Wextra -g -O2 -std=c++20 -DHITONANODE_LOCAL -DBENCHMARK main.cpp -o solver.out

format:
	find ./ -name "*.hpp" -o -name "*.cpp" | xargs clang-format --Werror -i -style=file

pyformat:
	uv run ruff format .
	uv run ruff check . --fix
	uv run mypy .

run:
	uv run -m heutopia.run

train:
	uv run -m heutopia.train

show-studies:
	uv run optuna studies  --storage sqlite:///train.db

clean:
	${RM} solver.out
