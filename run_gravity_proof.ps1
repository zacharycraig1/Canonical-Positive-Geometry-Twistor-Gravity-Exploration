param(
    [string]$LogFile = "RESULTS/gravity_proof_run.log"
)

# Create RESULTS directory if it doesn't exist
$resultsDir = "RESULTS"
if (-not (Test-Path $resultsDir)) {
    New-Item -ItemType Directory -Path $resultsDir | Out-Null
}

# Run the proof via Docker SageMath
Write-Host ("=" * 60)
Write-Host "Running 6-Point MHV Gravity Positive Geometry Proof"
Write-Host ("=" * 60)
Write-Host ""

$sageScript = "src/gravity_proof/run_proof.sage"

# Run via Docker and capture output
$command = "docker run --rm -v ${PWD}:/workspace sage-cursor sage $sageScript"

try {
    Write-Host "Executing: $command"
    Write-Host ""
    
    # Execute and capture both stdout and stderr
    $output = & docker run --rm -v ${PWD}:/workspace sage-cursor sage $sageScript 2>&1
    
    # Display output
    $output | ForEach-Object { Write-Host $_ }
    
    # Save to log file
    $output | Out-File -FilePath $LogFile -Encoding utf8
    
    Write-Host ""
    Write-Host ("=" * 60)
    Write-Host "Proof execution complete. Output saved to: $LogFile"
    Write-Host ("=" * 60)
    
    # Check if certificate was created
    $certFile = "RESULTS/gravity_proof_certificate.json"
    if (Test-Path $certFile) {
        Write-Host "Proof certificate: $certFile"
        
        # Try to read and display status
        try {
            $cert = Get-Content $certFile | ConvertFrom-Json
            Write-Host "Status: $($cert.status)"
            
            if ($cert.success_criteria) {
                Write-Host ""
                Write-Host "Success Criteria:"
                foreach ($key in $cert.success_criteria.PSObject.Properties.Name) {
                    $value = $cert.success_criteria.$key
                    $status = if ($value) { "PASS" } else { "FAIL" }
                    Write-Host "  $key : $status"
                }
            }
        } catch {
            Write-Host "Could not parse certificate JSON"
        }
    } else {
        Write-Host "Warning: Certificate file not found"
    }
    
    exit 0
    
} catch {
    Write-Host ""
    Write-Host "ERROR: Proof execution failed"
    Write-Host $_.Exception.Message
    Write-Host ""
    Write-Host "Full error output saved to: $LogFile"
    exit 1
}

