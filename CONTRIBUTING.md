# Contribution guidelines

## Package structure

Cacoa is the main class, providing the user interface. Most public methods inside belong to one of the two types: `estimate*` or `plot*`. The `plot*` methods don't include long computations, and these functions accepts only plotting parameters. The `estimate*` functions are the ones performing all real computations. Results of such computations are stored in the `test.results` field if they can be used for plotting right there, or in the `cache` if it's going to be used only in other `estimate*` functions. ***Note:*** it's still to be decided how do we deal with cache invalidation.

Cacoa operates with Conos, Seurat and possibly some other objects, which store experimental data. Such objects (one per Cacoa object) are stored in the `data.object` field. Extraction of the information from this object must be done through S3 `extract*` wrappers, stored in the "[access_wrappers.R](R/access_wrappers.R)" file. It's a good idea to keep `extract*` methods as short and fast as possible, but deviations are allowed.

Every `estimation*` function within Cacoa

## Code style

1. All functions are lowerCamelCase. See here why naming functions with dots is a bed idea (tldr: R uses dots for S3 functions)
2. All variables are lower case with dot as a separator (e.g. "n.pcs" or "count.matrix")
3. All files are named in a snake case with capital R as the extension (e.g. utility_functions.R)
4. Spaces around matrix operations, but not around function parameters (e.g. "x + 2 / 3" or "f(x=1, y=(2 / 3))")
5. Parentheses are required everywhere except one-line "return" or "stop" statements or short messages (i.e. "message" or "cat")
6. Trailing whitespaces must be removed before commiting files to git (see ["Strip trailing horizontal whitespace when saving" in RStudio](https://support.rstudio.com/hc/en-us/articles/200549016-Customizing-RStudio?mobile_site=true) and [Emacs](https://www.emacswiki.org/emacs/DeletingWhitespace#toc3) guides)
