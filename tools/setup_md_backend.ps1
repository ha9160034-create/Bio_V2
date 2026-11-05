# User-scoped, one-time WSL/MD setup. No Windows reboot or distribution removal.
# References: https://learn.microsoft.com/windows/wsl/basic-commands
param()
$ErrorActionPreference = 'Stop'
$taskStatus = Join-Path $PSScriptRoot 'md_backend_setup_status.json'
$taskLog = Join-Path $PSScriptRoot 'md_backend_setup.log'
$taskDistro = 'Ubuntu'

function Write-MDSetupState([string]$Stage, [string]$Message) {
    [ordered]@{
        stage = $Stage
        message = $Message
        distro = $taskDistro
        updated_at = [DateTime]::UtcNow.ToString('o')
        helper_pid = $PID
    } | ConvertTo-Json | Set-Content -LiteralPath $taskStatus -Encoding utf8
}
function Invoke-WSLSetup([string[]]$WslArguments) {
    & (Join-Path $env:SystemRoot 'System32\wsl.exe') @WslArguments 2>&1 |
        Out-File -LiteralPath $taskLog -Encoding utf8 -Append
    if ($LASTEXITCODE -ne 0) { throw "WSL command returned $LASTEXITCODE. See $taskLog." }
}

try {
    Write-MDSetupState 'checking_wsl' 'Checking the installed WSL runtime and distributions.'
    Invoke-WSLSetup @('--version')
    $taskInstalled = ((& wsl.exe --list --quiet 2>$null | Out-String) -replace "`0", '').Trim() -split '\r?\n'
    if ($taskDistro -notin $taskInstalled) {
        Write-MDSetupState 'installing_ubuntu' 'Downloading and registering Ubuntu without launching first-run prompts.'
        Invoke-WSLSetup @('--install', '-d', $taskDistro, '--no-launch', '--web-download')
    }
    Write-MDSetupState 'installing_md_tools' 'Installing GROMACS and ACPYPE/AmberTools inside Ubuntu.'
    $taskLinuxScript = '/mnt/c/Users/krar9/Desktop/Bio_V2-main/tools/install_md_backend.sh'
    Invoke-WSLSetup @('-d', $taskDistro, '-u', 'root', '--', 'bash', $taskLinuxScript)
    Write-MDSetupState 'backend_installed' 'Ubuntu and MD runtime packages installed. Verify project discovery and real MD next.'
} catch {
    Write-MDSetupState 'failed' $_.Exception.Message
    throw
}
