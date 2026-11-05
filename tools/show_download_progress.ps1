# Read-only progress display for the WSL installer downloaded by BioV2-WSL-Download.
$ErrorActionPreference = 'Stop'
$taskStatePath = Join-Path $PSScriptRoot 'wsl_setup_status.json'
$taskSetup = Get-Content -LiteralPath $taskStatePath -Raw | ConvertFrom-Json
$taskInstallerPath = $taskSetup.installer_path
$taskTotal = [long]$taskSetup.expected_bytes
if ($taskTotal -le 0) { throw 'The expected download size is unavailable.' }
$Host.UI.RawUI.WindowTitle = 'BioV2 - WSL download progress'
Write-Host 'WSL download progress (this window does not start another download).'
Write-Host 'Closing this window will not cancel the download.'
Write-Host ''
$taskPreviousBytes = $null
$taskPreviousTime = Get-Date
$taskSpeed = 0.0
while ($true) {
    $taskNow = Get-Date
    $taskBytes = if (Test-Path -LiteralPath $taskInstallerPath) {
        (Get-Item -LiteralPath $taskInstallerPath).Length
    } else { 0 }
    if ($null -ne $taskPreviousBytes) {
        $taskElapsed = ($taskNow - $taskPreviousTime).TotalSeconds
        if ($taskElapsed -gt 0) { $taskSpeed = [math]::Max(0, ($taskBytes - $taskPreviousBytes) / $taskElapsed / 1KB) }
    }
    $taskPercent = [math]::Min(100, 100.0 * $taskBytes / $taskTotal)
    $taskStatus = '{0:N1}% | {1:N1} / {2:N1} MB | {3:N0} KB/s' -f $taskPercent, ($taskBytes / 1MB), ($taskTotal / 1MB), $taskSpeed
    Write-Progress -Activity 'Downloading Microsoft WSL' -Status $taskStatus -PercentComplete ([int]$taskPercent)
    if ($taskBytes -ge $taskTotal) {
        Write-Progress -Activity 'Downloading Microsoft WSL' -Completed
        Write-Host 'Download complete. Installer verification and installation are separate steps.' -ForegroundColor Green
        break
    }
    $taskWorker = Get-ScheduledTask -TaskName 'BioV2-WSL-Download' -ErrorAction SilentlyContinue
    if ($null -eq $taskWorker -or $taskWorker.State -ne 'Running') {
        Write-Progress -Activity 'Downloading Microsoft WSL' -Completed
        Write-Host ('Download is not running. Current progress: ' + $taskStatus) -ForegroundColor Yellow
        break
    }
    $taskPreviousBytes = $taskBytes
    $taskPreviousTime = $taskNow
    Start-Sleep -Seconds 2
}
[void](Read-Host 'Press Enter to close')
