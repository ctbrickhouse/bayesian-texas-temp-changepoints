# 1. Render the Rmd to HTML and save the output path
html_output <- rmarkdown::render("analysis/report.Rmd")

# 2. Automatically open the file in your default browser
browseURL(html_output)


