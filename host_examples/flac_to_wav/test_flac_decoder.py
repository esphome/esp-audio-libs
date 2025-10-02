#!/usr/bin/env python3
"""
FLAC Decoder Test Script
Tests the esp-audio-libs FLAC decoder against ffmpeg for bit-perfect output
"""

import os
import sys
import subprocess
import hashlib
import json
from pathlib import Path
from datetime import datetime

# Configuration
FLAC_TO_WAV = "./flac_to_wav"
TEST_FILES_DIR = "flac-test-files"
RESULTS_DIR = "test_results"

# Test categories
TEST_CATEGORIES = ["subset", "uncommon", "faulty"]

class TestResult:
    def __init__(self, filename, category):
        self.filename = filename
        self.category = category
        self.our_decode_status = None
        self.our_decode_error = None
        self.ffmpeg_decode_status = None
        self.ffmpeg_decode_error = None
        self.comparison_result = None
        self.test_passed = None  # Primary result: based on MD5 verification
        self.pcm_match = None    # Secondary result: ffmpeg comparison
        self.our_md5 = None
        self.ffmpeg_md5 = None
        self.header_md5 = None
        self.md5_match = None

def run_command(cmd, timeout=30):
    """Run a command and return (success, stdout, stderr)"""
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            capture_output=True, 
            text=True, 
            timeout=timeout
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)

def extract_pcm_from_wav(wav_file):
    """Extract raw PCM data from WAV file (skip headers)"""
    try:
        with open(wav_file, 'rb') as f:
            # Read first 12 bytes to verify it's a WAV file
            riff = f.read(4)
            if riff != b'RIFF':
                return None
            
            f.read(4)  # File size
            wave = f.read(4)
            if wave != b'WAVE':
                return None
            
            # Find the data chunk
            while True:
                chunk_id = f.read(4)
                if not chunk_id:
                    return None
                    
                chunk_size = int.from_bytes(f.read(4), 'little')
                
                if chunk_id == b'data':
                    # Found data chunk, read PCM data
                    return f.read(chunk_size)
                else:
                    # Skip this chunk
                    f.seek(chunk_size, 1)
    except Exception as e:
        print(f"Error reading WAV file {wav_file}: {e}")
        return None

def compute_md5(data):
    """Compute MD5 hash of data"""
    if data is None:
        return None
    return hashlib.md5(data).hexdigest()

def get_flac_bit_depth(flac_file):
    """Get the bit depth of a FLAC file using ffprobe"""
    cmd = f'ffprobe -v error -select_streams a:0 -show_entries stream=bits_per_raw_sample -of default=noprint_wrappers=1:nokey=1 "{flac_file}"'
    success, stdout, stderr = run_command(cmd)
    if success and stdout.strip():
        try:
            return int(stdout.strip())
        except:
            pass
    return None

def get_ffmpeg_codec_for_bit_depth(bit_depth):
    """Get the appropriate ffmpeg codec for a given bit depth"""
    if bit_depth is None:
        return ""  # Let ffmpeg decide
    elif bit_depth <= 8:
        return "-c:a pcm_u8"
    elif bit_depth <= 16:
        return "-c:a pcm_s16le"
    elif bit_depth <= 24:
        return "-c:a pcm_s24le"
    elif bit_depth <= 32:
        return "-c:a pcm_s32le"
    else:
        return ""  # Let ffmpeg decide

def extract_md5_from_output(stdout):
    """Extract MD5 signature from decoder output"""
    for line in stdout.splitlines():
        if "MD5 signature:" in line:
            # Extract the hex string after "MD5 signature:"
            md5_str = line.split("MD5 signature:")[1].strip()
            return md5_str
    return None

def extract_md5_verification_from_output(stdout):
    """Extract MD5 verification result from decoder output
    Returns (expected_md5, computed_md5, passed) or (None, None, None) if not found
    """
    expected = None
    computed = None
    passed = None

    lines = stdout.splitlines()
    for i, line in enumerate(lines):
        if "Expected MD5:" in line:
            expected = line.split("Expected MD5:")[1].strip()
        elif "Computed MD5:" in line:
            computed = line.split("Computed MD5:")[1].strip()
        elif "Result: PASS" in line:
            passed = True
        elif "Result: FAIL" in line:
            passed = False
        elif "Status: SKIPPED" in line:
            # MD5 not available in file
            return None, None, None

    return expected, computed, passed

def test_single_file(flac_file, category):
    """Test a single FLAC file"""
    result = TestResult(os.path.basename(flac_file), category)

    # Get bit depth to use appropriate ffmpeg codec
    bit_depth = get_flac_bit_depth(flac_file)

    # Create output paths
    base_name = os.path.splitext(os.path.basename(flac_file))[0]
    our_output = os.path.join(RESULTS_DIR, category, "our_decoder", f"{base_name}.wav")
    ffmpeg_output = os.path.join(RESULTS_DIR, category, "ffmpeg", f"{base_name}.wav")

    # Ensure output directories exist
    os.makedirs(os.path.dirname(our_output), exist_ok=True)
    os.makedirs(os.path.dirname(ffmpeg_output), exist_ok=True)
    
    # Test our decoder
    cmd = f'"{FLAC_TO_WAV}" "{flac_file}" "{our_output}"'
    success, stdout, stderr = run_command(cmd)
    result.our_decode_status = "success" if success else "failed"
    if not success:
        result.our_decode_error = stderr.strip() if stderr else stdout.strip()
    else:
        # Extract MD5 verification from decoder output
        expected_md5, computed_md5, md5_passed = extract_md5_verification_from_output(stdout)
        result.header_md5 = expected_md5
        result.our_md5 = computed_md5
        if md5_passed is not None:
            result.md5_match = md5_passed
    
    # Test ffmpeg decoder with appropriate codec for bit depth
    codec_arg = get_ffmpeg_codec_for_bit_depth(bit_depth)
    cmd = f'ffmpeg -i "{flac_file}" {codec_arg} -f wav -y "{ffmpeg_output}" 2>&1'
    success, stdout, stderr = run_command(cmd)
    result.ffmpeg_decode_status = "success" if success else "failed"
    if not success:
        result.ffmpeg_decode_error = (stderr + stdout).strip()
    
    # Determine test result based on MD5 verification (primary)
    if result.our_decode_status == "success":
        if result.md5_match is True:
            # MD5 verification passed - this is definitive proof of correct decode
            result.test_passed = True
            result.comparison_result = "PASS - MD5 verified"
        elif result.md5_match is False:
            # MD5 verification failed - decode is incorrect
            result.test_passed = False
            result.comparison_result = "FAIL - MD5 mismatch"
        elif result.md5_match is None:
            # No MD5 in file - fall back to ffmpeg comparison
            if result.ffmpeg_decode_status == "success":
                our_pcm = extract_pcm_from_wav(our_output)
                ffmpeg_pcm = extract_pcm_from_wav(ffmpeg_output)

                if our_pcm is not None and ffmpeg_pcm is not None:
                    result.ffmpeg_md5 = compute_md5(ffmpeg_pcm)
                    if our_pcm == ffmpeg_pcm:
                        result.pcm_match = True
                        result.test_passed = True
                        result.comparison_result = "PASS - Matches ffmpeg (no MD5 in file)"
                    else:
                        result.pcm_match = False
                        result.test_passed = False
                        result.comparison_result = f"FAIL - PCM mismatch with ffmpeg (our: {len(our_pcm)} bytes, ffmpeg: {len(ffmpeg_pcm)} bytes)"
                else:
                    result.test_passed = None
                    result.comparison_result = "ERROR - Could not extract PCM data"
            else:
                # Our decoder succeeded but no MD5 and ffmpeg failed
                result.test_passed = None
                result.comparison_result = "UNKNOWN - No MD5 available and ffmpeg failed"

        # Secondary check: compare with ffmpeg if both succeeded (for additional validation)
        if result.ffmpeg_decode_status == "success" and result.md5_match is not None:
            our_pcm = extract_pcm_from_wav(our_output)
            ffmpeg_pcm = extract_pcm_from_wav(ffmpeg_output)

            if our_pcm is not None and ffmpeg_pcm is not None:
                result.ffmpeg_md5 = compute_md5(ffmpeg_pcm)
                result.pcm_match = (our_pcm == ffmpeg_pcm)

                # Add note if ffmpeg disagrees with MD5 result
                if result.test_passed and not result.pcm_match:
                    result.comparison_result += " (WARNING: ffmpeg output differs)"
                elif not result.test_passed and result.pcm_match:
                    result.comparison_result += " (NOTE: ffmpeg output matches despite MD5 fail)"
                elif result.test_passed and result.pcm_match:
                    result.comparison_result += " + matches ffmpeg"

    elif result.our_decode_status == "failed":
        result.test_passed = False
        if result.ffmpeg_decode_status == "failed":
            result.comparison_result = "EXPECTED - Both decoders failed (likely invalid file)"
            result.test_passed = None  # Expected failure, not a real failure
        else:
            result.comparison_result = "FAIL - Decoder failed but ffmpeg succeeded"
    else:
        result.test_passed = None
        result.comparison_result = "ERROR - Decode status unknown"
    
    # Clean up output files to save space (optional)
    # Uncomment if you want to delete the WAV files after comparison
    # if os.path.exists(our_output):
    #     os.remove(our_output)
    # if os.path.exists(ffmpeg_output):
    #     os.remove(ffmpeg_output)
    
    return result

def test_category(category):
    """Test all files in a category"""
    category_dir = os.path.join(TEST_FILES_DIR, category)
    results = []
    
    # Find all FLAC files in the category
    flac_files = sorted(Path(category_dir).glob("*.flac"))
    
    print(f"\nTesting {category} files ({len(flac_files)} files)...")
    
    for i, flac_file in enumerate(flac_files, 1):
        print(f"  [{i}/{len(flac_files)}] Testing {flac_file.name}...", end="")
        result = test_single_file(str(flac_file), category)
        results.append(result)

        # Print immediate result based on test_passed
        if result.test_passed is True:
            print(" ✓ PASS")
        elif result.test_passed is False:
            print(" ✗ FAIL")
        else:
            print(f" - {result.comparison_result[:40]}")
    
    return results

def generate_report(all_results):
    """Generate test report"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Calculate statistics
    stats = {
        "total": len(all_results),
        "passed": sum(1 for r in all_results if r.test_passed is True),
        "failed": sum(1 for r in all_results if r.test_passed is False),
        "errors": sum(1 for r in all_results if r.test_passed is None),
    }
    
    # Generate text report
    report_lines = [
        "=" * 80,
        "FLAC Decoder Test Report",
        f"Generated: {timestamp}",
        "=" * 80,
        "",
        "SUMMARY",
        "-" * 40,
        f"Total files tested: {stats['total']}",
        f"Passed (bit-perfect): {stats['passed']} ({stats['passed']*100//stats['total'] if stats['total'] else 0}%)",
        f"Failed: {stats['failed']}",
        f"Errors/Expected failures: {stats['errors']}",
        "",
    ]
    
    # Detailed results by category
    for category in TEST_CATEGORIES:
        category_results = [r for r in all_results if r.category == category]
        if not category_results:
            continue
            
        report_lines.extend([
            "",
            f"{category.upper()} FILES ({len(category_results)} files)",
            "-" * 40,
        ])
        
        for result in category_results:
            status = "✓" if result.test_passed is True else "✗" if result.test_passed is False else "?"
            report_lines.append(f"{status} {result.filename}: {result.comparison_result}")

            # Show MD5 details if there was a mismatch
            if result.md5_match is False and result.header_md5:
                report_lines.append(f"    Expected MD5: {result.header_md5}")
                report_lines.append(f"    Computed MD5: {result.our_md5}")

            # Show ffmpeg comparison details if relevant
            if result.pcm_match is not None:
                ffmpeg_status = "matches" if result.pcm_match else "differs"
                report_lines.append(f"    ffmpeg comparison: {ffmpeg_status}")

            if result.our_decode_error:
                report_lines.append(f"    Our decoder error: {result.our_decode_error}")
            if result.ffmpeg_decode_error and "ffmpeg failed" in result.comparison_result:
                report_lines.append(f"    ffmpeg error: {result.ffmpeg_decode_error[:100]}")
    
    # Ensure results directory exists
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Write text report
    report_text = "\n".join(report_lines)
    report_file = os.path.join(RESULTS_DIR, "test_report.txt")
    with open(report_file, "w") as f:
        f.write(report_text)
    
    # Also save JSON report for programmatic access
    json_report = {
        "timestamp": timestamp,
        "summary": stats,
        "results": [
            {
                "file": r.filename,
                "category": r.category,
                "our_status": r.our_decode_status,
                "ffmpeg_status": r.ffmpeg_decode_status,
                "test_passed": r.test_passed,
                "pcm_match": r.pcm_match,
                "md5_match": r.md5_match,
                "result": r.comparison_result,
                "header_md5": r.header_md5,
                "our_md5": r.our_md5,
                "ffmpeg_md5": r.ffmpeg_md5,
            }
            for r in all_results
        ]
    }
    
    json_file = os.path.join(RESULTS_DIR, "test_report.json")
    with open(json_file, "w") as f:
        json.dump(json_report, f, indent=2)
    
    return report_text, report_file, json_file

def main():
    """Main test runner"""
    print("FLAC Decoder Test Suite")
    print("=" * 40)
    
    # Check prerequisites
    if not os.path.exists(FLAC_TO_WAV):
        print(f"Error: flac_to_wav not found at {FLAC_TO_WAV}")
        print("Please build it first in host_examples/flac_to_wav/")
        return 1
    
    if not os.path.exists(TEST_FILES_DIR):
        print(f"Error: Test files not found at {TEST_FILES_DIR}")
        print("Please clone https://github.com/ietf-wg-cellar/flac-test-files")
        return 1
    
    # Check ffmpeg
    success, _, _ = run_command("ffmpeg -version")
    if not success:
        print("Error: ffmpeg not found. Please install ffmpeg.")
        return 1
    
    print(f"✓ Found flac_to_wav at {FLAC_TO_WAV}")
    print(f"✓ Found test files at {TEST_FILES_DIR}")
    print(f"✓ ffmpeg is available")
    print()
    
    # Run tests for each category
    all_results = []
    for category in TEST_CATEGORIES:
        results = test_category(category)
        all_results.extend(results)
    
    # Generate report
    print("\nGenerating report...")
    report_text, report_file, json_file = generate_report(all_results)
    
    # Print summary
    print("\n" + "=" * 40)
    print("TEST COMPLETE")
    print("=" * 40)
    
    stats = {
        "total": len(all_results),
        "passed": sum(1 for r in all_results if r.test_passed is True),
        "failed": sum(1 for r in all_results if r.test_passed is False),
        "errors": sum(1 for r in all_results if r.test_passed is None),
    }

    print(f"Total: {stats['total']} files")
    print(f"Passed: {stats['passed']} ({stats['passed']*100//stats['total'] if stats['total'] else 0}%)")
    print(f"Failed: {stats['failed']}")
    print(f"Errors/Expected: {stats['errors']}")
    print()
    print(f"Full report saved to: {report_file}")
    print(f"JSON report saved to: {json_file}")

    # Return non-zero if any tests failed
    return 0 if stats['failed'] == 0 else 1

if __name__ == "__main__":
    sys.exit(main())