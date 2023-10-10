#!/usr/bin/env pwsh

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$rootDir = Join-Path -Path $scriptDir -ChildPath ..\..\

& "$rootDir/bootstrap-conda"

$platforms = @{
    "linux-64" = "linux"
    "linux-aarch64" = "linux-aarch"
    "osx-64" = "macos"
    "osx-arm64" = "macos-arm"
}
$pythons = "3.9", "3.10", "3.11"
$tags = "", "-dev"
$sources = "", "src"

foreach ($platform in $platforms.GetEnumerator()) {
    foreach ($python in $pythons) {
        foreach ($tag in $tags) {
            foreach ($src in $sources) {
                $envFile = Join-Path -Path $rootDir -ChildPath "$src\environment$tag-$python.yml"
                $lockFile = Join-Path -Path $rootDir -ChildPath "$src\environment$tag-$python-$($platform.Value)"

                if (-not (Test-Path $envFile)) {
                    continue
                }

                echo "Updating lock file for $envFile at $lockFile"
                & "conda-lock" --channel conda-forge --kind env --platform $platform.Key --file $envFile --lockfile $lockFile --filename-template $lockFile
            }
        }
    }
}
