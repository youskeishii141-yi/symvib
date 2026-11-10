# GUI 版の配布用 exe をビルドする (Windows PowerShell)。
#   使い方:  .\build_exe.ps1
#   出力:    dist\SymVib\  (フォルダごと zip にして配布)
#            dist\SymVib.zip
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

if (-not (Test-Path ".venv")) {
    python -m venv .venv
}
$py = ".\.venv\Scripts\python.exe"
& $py -m pip install --upgrade pip
& $py -m pip install -r requirements.txt -r requirements-dev.txt

# ビルド前にテストを通す
& $py -m pytest -q tests
if ($LASTEXITCODE -ne 0) { throw "テストが失敗したためビルドを中止しました" }

& $py -m PyInstaller SymVib.spec --noconfirm --clean
if ($LASTEXITCODE -ne 0) { throw "PyInstaller が失敗しました" }

# 同梱ライブラリのライセンス表示 (再配布の条件) と本体のライセンス・説明書を添える
& $py tools\third_party_notices.py dist\SymVib\THIRD_PARTY_LICENSES.txt
if ($LASTEXITCODE -ne 0) { throw "THIRD_PARTY_LICENSES.txt の作成に失敗しました" }
Copy-Item LICENSE, README.md, README.ja.md -Destination dist\SymVib\

# 起動時の自己診断に使うファイル一覧を記録し、exe の中で解析まで動くか確かめる
& $py tools\bundle_manifest.py dist\SymVib
$out = Join-Path $env:TEMP "symvib_self_test.txt"
$p = Start-Process -FilePath "dist\SymVib\SymVib.exe" -ArgumentList "--self-test", $out -PassThru -Wait
Get-Content $out
if ($p.ExitCode -ne 0) { throw "exe の self-test が失敗しました" }

Compress-Archive -Path "dist\SymVib" -DestinationPath "dist\SymVib.zip" -Force
Write-Host "完成: dist\SymVib\SymVib.exe"
Write-Host "配布用: dist\SymVib.zip"
