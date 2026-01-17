jupyter nbconvert --to markdown .ipynb
jupyter nbconvert --to markdown plot_SEC.ipynb

$env:HTTP_PROXY = "http://127.0.0.1:7890"
$env:HTTPS_PROXY = "http://127.0.0.1:7890"

git add .
git commit -m ""
git push origin main