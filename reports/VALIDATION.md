# Website validation

Validated on 8 September 2026 against the executed Python notebooks in the accompanying content checkout.

- Jekyll 4.3.3 builds the complete site successfully with Ruby 3.1.3.
- All 58 distinct Python tours have a reading page, an actual computed figure preview, an executed notebook download, and a Colab link.
- The exporter requires passing execution reports matching the exact source notebook hashes; the release export did not use the development `--allow-stale` option.
- The actual search script passes checks for queries, multiple words, accents, topic filters, empty results, reset, URL state, and literal input.
- All primary navigation routes and all local links and image references in the 58 reading pages resolve in the built site.
- Downloaded notebooks contain no execution errors and every code cell has an execution count.
- Both repository diffs pass whitespace checks.

The content validation report records **62/62 notebooks**, **1,892 executed code cells**, and **71 passing tests**. Four historical variants are retained in the content branch and omitted from the website catalogue. The unfinished `python/todo/` conversions and legacy-language runtimes were not execution-tested.

Browser visual testing and remote Colab execution were not performed. Colab buttons reference `master`; publish the content changes before publishing this website branch. GitHub Pages serves the root of `gh-pages` at `www.numerical-tours.com`.

Reproduce the checks with the instructions in [README.md](../README.md).
