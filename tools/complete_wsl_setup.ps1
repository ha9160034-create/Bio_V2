#Requires -RunAsAdministrator
<#
Finish the user-authorized WSL installation after the official x64 MSI downloads.
References: https://learn.microsoft.com/windows/wsl/install
            https://learn.microsoft.com/windows/wsl/install-manual
Never restarts Windows, removes distributions, or changes security policies.
#>
param(
    [string]$Installer = 'C:\Users\krar9\AppData\Local\Temp\BioV2-WSL-6ce02f48bb4c41ec8de32f124cb5ac8c\wsl.2.7.14.0.x64.msi',
    [int]$BundleInstallerPid = 0
)
$ErrorActionPreference = 'Stop'
$taskStatePath = Join-Path $PSScriptRoot 'wsl_setup_status.json'
$taskTranscript = Join-Path $PSScriptRoot 'complete_wsl_setup.log'
$taskExpectedSize = 258990080
$taskExpectedHash = 'DB084E536279A59E90A26EC598D8AA8A4DFF8309F41D078FD06242953AC1EBCD'

function Write-SetupState([string]$Stage, [string]$Message) {
    $taskBytes = if (Test-Path -LiteralPath $Installer) { (Get-Item -LiteralPath $Installer).Length } else { 0 }
    [ordered]@{
        stage = $Stage
        message = $Message
        updated_at = [DateTime]::UtcNow.ToString('o')
        installer_path = $Installer
        downloaded_bytes = $taskBytes
        expected_bytes = $taskExpectedSize
        helper_pid = $PID
    } | ConvertTo-Json | Set-Content -LiteralPath $taskStatePath -Encoding utf8
    Write-Output "$Stage : $Message"
}

Start-Transcript -Path $taskTranscript -Force
try {
    # Only cancel the exact redundant WSL downloader started by this chat.
    $taskBundle = if ($BundleInstallerPid -gt 0) {
        Get-CimInstance Win32_Process -Filter "ProcessId = $BundleInstallerPid"
    } else { $null }
    if ($taskBundle) {
        if ($taskBundle.Name -ne 'winget.exe' -or $taskBundle.CommandLine -notlike '*--id Microsoft.WSL*') {
            throw 'The redundant installer PID no longer belongs to this WSL installation.'
        }
        Stop-Process -Id $BundleInstallerPid -ErrorAction Stop
    }

    $taskDeadline = [DateTime]::UtcNow.AddHours(8)
    while ($true) {
        if ([DateTime]::UtcNow -gt $taskDeadline) { throw 'The WSL download exceeded eight hours.' }
        $taskBytes = if (Test-Path -LiteralPath $Installer) { (Get-Item -LiteralPath $Installer).Length } else { 0 }
        if ($taskBytes -gt $taskExpectedSize) { throw 'Unexpected WSL installer size.' }
        if ($taskBytes -eq $taskExpectedSize) { break }
        Write-SetupState 'downloading' 'Waiting for the official Microsoft x64 MSI download.'
        Start-Sleep -Seconds 15
    }

    Write-SetupState 'verifying' 'Verifying the release SHA256 and Microsoft digital signature.'
    $taskHash = (Get-FileHash -LiteralPath $Installer -Algorithm SHA256).Hash
    if ($taskHash -ne $taskExpectedHash) { throw 'WSL installer SHA256 does not match the Microsoft release.' }
    $taskSignature = Get-AuthenticodeSignature -LiteralPath $Installer
    if ($taskSignature.Status -ne 'Valid' -or $taskSignature.SignerCertificate.Subject -notmatch 'O=Microsoft Corporation') {
        throw 'WSL installer does not have a valid Microsoft signature.'
    }

    Write-SetupState 'installing_wsl' 'Installing the verified WSL runtime without restarting Windows.'
    $taskMsiLog = Join-Path $PSScriptRoot 'wsl_msi_install.log'
    $taskMsiArgs = "/i `"$Installer`" /qn /norestart /L*v `"$taskMsiLog`""
    $taskMsi = Start-Process -FilePath (Join-Path $env:SystemRoot 'System32\msiexec.exe') -ArgumentList $taskMsiArgs -Wait -PassThru -WindowStyle Hidden
    if ($taskMsi.ExitCode -notin @(0, 3010)) { throw "WSL MSI installation failed: $($taskMsi.ExitCode)." }

    $taskRestart = $taskMsi.ExitCode -eq 3010
    foreach ($taskFeatureName in @('Microsoft-Windows-Subsystem-Linux', 'VirtualMachinePlatform')) {
        $taskFeature = Get-WindowsOptionalFeature -Online -FeatureName $taskFeatureName
        if ($taskFeature.State -ne 'Enabled') {
            Write-SetupState 'enabling_features' "Enabling the WSL prerequisite $taskFeatureName without restarting."
            $taskEnabled = Enable-WindowsOptionalFeature -Online -FeatureName $taskFeatureName -All -NoRestart
            $taskRestart = $taskRestart -or $taskEnabled.RestartNeeded
        }
    }
    if ($taskRestart) {
        Write-SetupState 'restart_required' 'WSL is installed. Save your work and restart Windows before continuing Ubuntu/GROMACS setup.'
        exit 0
    }

    & wsl.exe --version
    if ($LASTEXITCODE -ne 0) { throw 'WSL version check failed after installing its runtime.' }
    Write-SetupState 'installing_ubuntu' 'Installing Ubuntu without first-run user prompts or an automatic restart.'
    & wsl.exe --install -d Ubuntu --no-launch --web-download
    if ($LASTEXITCODE -ne 0) { throw "Ubuntu setup returned exit code $LASTEXITCODE; inspect the transcript." }
    Write-SetupState 'ubuntu_installed' 'WSL and Ubuntu installation finished. Verify the distribution and install the MD backend next.'
} catch {
    Write-SetupState 'failed' $_.Exception.Message
    throw
} finally {
    Stop-Transcript
}
