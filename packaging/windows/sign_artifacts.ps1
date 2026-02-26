param(
    [Parameter(Mandatory = $true)]
    [string]$InputDir,
    [string]$TimestampUrl = "http://timestamp.digicert.com"
)

$certBase64 = $env:WINDOWS_CODESIGN_CERT_BASE64
$certPassword = $env:WINDOWS_CODESIGN_CERT_PASSWORD

if ([string]::IsNullOrWhiteSpace($certBase64) -or [string]::IsNullOrWhiteSpace($certPassword)) {
    Write-Host "No se configuraron secretos de firma Authenticode; se omite firmado."
    exit 0
}

if (-not (Test-Path -Path $InputDir)) {
    throw "No existe InputDir: $InputDir"
}

$targets = Get-ChildItem -Path $InputDir -File | Where-Object { $_.Extension -in @(".exe", ".msi") }
if ($null -eq $targets -or $targets.Count -eq 0) {
    Write-Host "No hay artefactos .exe/.msi para firmar en $InputDir."
    exit 0
}

$signtool = (Get-Command signtool.exe -ErrorAction Stop).Source
$tempPfx = Join-Path ($env:RUNNER_TEMP ? $env:RUNNER_TEMP : $env:TEMP) "chemuson_codesign.pfx"

try {
    [System.IO.File]::WriteAllBytes($tempPfx, [Convert]::FromBase64String($certBase64))
    foreach ($file in $targets) {
        Write-Host "Firmando $($file.Name)"
        & $signtool sign `
            /f $tempPfx `
            /p $certPassword `
            /fd SHA256 `
            /td SHA256 `
            /tr $TimestampUrl `
            "$($file.FullName)"
        if ($LASTEXITCODE -ne 0) {
            throw "Fallo de firma en $($file.FullName)"
        }
    }
}
finally {
    if (Test-Path -Path $tempPfx) {
        Remove-Item -Path $tempPfx -Force
    }
}
