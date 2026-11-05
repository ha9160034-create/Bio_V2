# Read-only display of the independent Ubuntu / MD installation worker.
$ErrorActionPreference = 'Stop'
$Host.UI.RawUI.WindowTitle = 'BioV2 - Ubuntu and MD installation'
$taskStatusPath = Join-Path $PSScriptRoot 'md_backend_setup_status.json'
$taskLogPath = Join-Path $PSScriptRoot 'md_backend_setup.log'
$taskLastStage = ''
$taskUbuntuFile = Join-Path ([IO.Path]::GetTempPath()) 'Ubuntu.wsl'
# Size obtained from the official WSL distribution manifest and Ubuntu HTTP headers.
$taskUbuntuTotal = 418495746
$taskLastBytes = $null
$taskLastTime = Get-Date
Clear-Host
Write-Host ('=' * 68) -ForegroundColor DarkCyan
Write-Host '        Bio_V2  |  Ubuntu + MD backend installation' -ForegroundColor Cyan
Write-Host ('=' * 68) -ForegroundColor DarkCyan
Write-Host '  هذه النافذة للعرض فقط؛ إغلاقها لا يوقف التثبيت.' -ForegroundColor Gray
Write-Host '  المرحلة الحالية والتقدم يظهران بالأسفل.' -ForegroundColor Gray
Write-Host ''
Write-Host ''
while ($true) {
    try {
        $taskStatus = Get-Content -LiteralPath $taskStatusPath -Raw | ConvertFrom-Json
    } catch {
        Start-Sleep -Seconds 2
        continue
    }
    if ($taskStatus.stage -ne $taskLastStage) {
        Write-Host ''
        Write-Host ('-' * 68) -ForegroundColor DarkGray
        Write-Host ("  [{0}]  {1}" -f (Get-Date -Format 'HH:mm:ss'), $taskStatus.message) -ForegroundColor Cyan
        $taskLastStage = $taskStatus.stage
    }
    $taskWorker = Get-ScheduledTask -TaskName 'BioV2-MD-Backend-Setup' -ErrorAction SilentlyContinue
    if ($taskStatus.stage -eq 'failed') {
        Write-Progress -Activity 'Ubuntu and MD setup' -Completed
        Write-Host $taskStatus.message -ForegroundColor Red
        break
    }
    if ($taskStatus.stage -eq 'backend_installed') {
        Write-Progress -Activity 'Ubuntu and MD setup' -Completed
        Write-Host 'Packages installed. Real MD verification is a separate next step.' -ForegroundColor Green
        break
    }
    if ($null -eq $taskWorker -or $taskWorker.State -ne 'Running') {
        Write-Progress -Activity 'Ubuntu and MD setup' -Completed
        Write-Host 'The setup worker is not running. Check the setup log.' -ForegroundColor Yellow
        break
    }
    if ($taskStatus.stage -eq 'installing_ubuntu' -and (Test-Path -LiteralPath $taskUbuntuFile)) {
        $taskNow = Get-Date
        $taskBytes = (Get-Item -LiteralPath $taskUbuntuFile).Length
        $taskRate = 0.0
        $taskElapsed = ($taskNow - $taskLastTime).TotalSeconds
        if ($null -ne $taskLastBytes -and $taskElapsed -gt 0) {
            $taskRate = [math]::Max(0, ($taskBytes - $taskLastBytes) / $taskElapsed / 1000)
        }
        $taskPercent = [math]::Min(100, 100.0 * $taskBytes / $taskUbuntuTotal)
        $taskProgress = '{0:N1}% | {1:N1} / {2:N1} MB | {3:N0} KB/s' -f $taskPercent, ($taskBytes / 1000000), ($taskUbuntuTotal / 1000000), $taskRate
        Write-Progress -Activity 'Downloading Ubuntu' -Status $taskProgress -PercentComplete ([int]$taskPercent)
        $taskLastBytes = $taskBytes
        $taskLastTime = $taskNow
    } else {
        $taskTail = if (Test-Path -LiteralPath $taskLogPath) {
            ((Get-Content -LiteralPath $taskLogPath -Tail 8 | ForEach-Object { $_ -replace "`0", '' } |
                Where-Object { $_ -and $_ -notmatch 'version:|WSLg|MSRDC|Direct3D|DXCore' }) -join ' ')
        } else { '' }
        Write-Progress -Activity 'Ubuntu and MD setup is running' -Status $taskStatus.message -CurrentOperation $taskTail
        if ($taskTail) {
            Write-Host ("  آخر نشاط: {0}" -f $taskTail) -ForegroundColor DarkGray
        }
    }
    Start-Sleep -Seconds 2
}
if (Test-Path -LiteralPath $taskLogPath) {
    Write-Host ''
    Write-Host (((Get-Content -LiteralPath $taskLogPath -Tail 12) -join "`n") -replace "`0", '')
}
[void](Read-Host 'Press Enter to close')
