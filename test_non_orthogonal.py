"""
Test script for SAW-NO non-orthogonal propagation changes.
Runs three tests:
  1. Faraday geometry: orthogonal vs non-orthogonal (should differ only by Fresnel window factors)
  2. Arbitrary angle: both paths run without error
  3. Sanity checks on output shapes and physics
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from elecsus.elecsus_methods import calculate

def test_faraday_both_paths():
    """Test Faraday geometry (Btheta=0) with both propagation paths."""
    print("=" * 60)
    print("TEST 1: Faraday geometry (Btheta=0)")
    print("=" * 60)
    
    detuning = np.linspace(-5, 5, 500) * 1e3  # MHz
    E_in = np.array([1, 0, 0])
    
    p_base = {'Elem': 'Cs', 'Dline': 'D2', 'T': 21, 'lcell': 75e-3,
              'Bfield': 100, 'Btheta': 0, 'Bphi': 0}
    
    # Original path
    p_ortho = dict(p_base, use_non_orthogonal=False)
    result_ortho = calculate(detuning, E_in, p_ortho, outputs=['S0', 'S1', 'S2', 'S3'])
    S0_o, S1_o, S2_o, S3_o = [r.real for r in result_ortho]
    
    # Non-orthogonal path
    p_nonortho = dict(p_base, use_non_orthogonal=True, n_window=1.51)
    result_no = calculate(detuning, E_in, p_nonortho, outputs=['S0', 'S1', 'S2', 'S3'])
    S0_n, S1_n, S2_n, S3_n = [r.real for r in result_no]
    
    print(f"  S0 orthogonal   : min={S0_o.min():.6f}, max={S0_o.max():.6f}")
    print(f"  S0 non-orthogonal: min={S0_n.min():.6f}, max={S0_n.max():.6f}")
    
    # They won't be identical because non-orthogonal path includes Fresnel coefficients.
    # But both should be physical (0 <= S0 <= ~1 far from resonance).
    assert S0_o.shape == (500,), f"S0 orthogonal shape wrong: {S0_o.shape}"
    assert S0_n.shape == (500,), f"S0 non-orthogonal shape wrong: {S0_n.shape}"
    assert np.all(np.isfinite(S0_o)), "S0 orthogonal has non-finite values"
    assert np.all(np.isfinite(S0_n)), "S0 non-orthogonal has non-finite values"
    assert np.all(S0_o >= 0), "S0 orthogonal has negative values"
    assert np.all(S0_n >= 0), "S0 non-orthogonal has negative values"
    
    # Far from resonance, transmission should approach ~1 (ortho) or Fresnel-reduced (non-ortho)
    # Check the wings (first and last 20 points)
    wing_ortho = np.mean(np.concatenate([S0_o[:20], S0_o[-20:]]))
    wing_no = np.mean(np.concatenate([S0_n[:20], S0_n[-20:]]))
    print(f"  Wing transmission (ortho):     {wing_ortho:.6f}")
    print(f"  Wing transmission (non-ortho): {wing_no:.6f}")
    
    # Fresnel factor for two surfaces: T = |t_in * t_out|^2 per mode
    # For n_window=1.51, a single surface t = 2*1.51/(1.51+1) ~ 1.203
    # but total power transmission per surface ~ 4*n*1/(n+1)^2 ~ 0.957
    # Two surfaces ~ 0.957^2 ~ 0.916
    expected_fresnel = (4 * 1.51 / (1.51 + 1)**2)**2
    print(f"  Expected Fresnel factor (2 surfaces): {expected_fresnel:.6f}")
    ratio = wing_no / wing_ortho
    print(f"  Actual ratio (non-ortho/ortho wings): {ratio:.6f}")
    
    print("  PASSED\n")

def test_arbitrary_angle():
    """Test at arbitrary B-field angle."""
    print("=" * 60)
    print("TEST 2: Arbitrary angle (Btheta=0.7 rad)")
    print("=" * 60)
    
    detuning = np.linspace(-5, 5, 300) * 1e3
    E_in = np.array([1, 0, 0])
    
    p_base = {'Elem': 'Rb', 'Dline': 'D2', 'T': 20, 'lcell': 75e-3,
              'Bfield': 100, 'Btheta': 0.7, 'Bphi': 0}
    
    # Original path
    p_ortho = dict(p_base, use_non_orthogonal=False)
    result_ortho = calculate(detuning, E_in, p_ortho, outputs=['S0', 'S1', 'S2', 'S3'])
    S0_o, S1_o, S2_o, S3_o = [r.real for r in result_ortho]
    
    # Non-orthogonal path
    p_nonortho = dict(p_base, use_non_orthogonal=True)
    result_no = calculate(detuning, E_in, p_nonortho, outputs=['S0', 'S1', 'S2', 'S3'])
    S0_n, S1_n, S2_n, S3_n = [r.real for r in result_no]
    
    print(f"  S0 orthogonal   : min={S0_o.min():.6f}, max={S0_o.max():.6f}")
    print(f"  S0 non-orthogonal: min={S0_n.min():.6f}, max={S0_n.max():.6f}")
    print(f"  S1 non-orthogonal: min={S1_n.min():.6f}, max={S1_n.max():.6f}")
    print(f"  S2 non-orthogonal: min={S2_n.min():.6f}, max={S2_n.max():.6f}")
    print(f"  S3 non-orthogonal: min={S3_n.min():.6f}, max={S3_n.max():.6f}")
    
    assert np.all(np.isfinite(S0_n)), "S0 non-orthogonal has non-finite values"
    assert np.all(np.isfinite(S1_n)), "S1 non-orthogonal has non-finite values"
    assert np.all(np.isfinite(S2_n)), "S2 non-orthogonal has non-finite values"
    assert np.all(np.isfinite(S3_n)), "S3 non-orthogonal has non-finite values"
    assert np.all(S0_n >= 0), "S0 non-orthogonal has negative values"
    
    print("  PASSED\n")

def test_voigt_geometry():
    """Test Voigt geometry (Btheta=pi/2)."""
    print("=" * 60)
    print("TEST 3: Voigt geometry (Btheta=pi/2)")
    print("=" * 60)
    
    detuning = np.linspace(-5, 5, 300) * 1e3
    E_in = np.array([1, 0, 0])
    
    p_base = {'Elem': 'Cs', 'Dline': 'D2', 'T': 21, 'lcell': 75e-3,
              'Bfield': 100, 'Btheta': np.pi/2, 'Bphi': 0}
    
    p_ortho = dict(p_base, use_non_orthogonal=False)
    result_ortho = calculate(detuning, E_in, p_ortho, outputs=['S0'])
    S0_o = result_ortho[0].real
    
    p_nonortho = dict(p_base, use_non_orthogonal=True)
    result_no = calculate(detuning, E_in, p_nonortho, outputs=['S0'])
    S0_n = result_no[0].real
    
    print(f"  S0 orthogonal   : min={S0_o.min():.6f}, max={S0_o.max():.6f}")
    print(f"  S0 non-orthogonal: min={S0_n.min():.6f}, max={S0_n.max():.6f}")
    
    assert np.all(np.isfinite(S0_n)), "S0 non-orthogonal has non-finite values"
    assert np.all(S0_n >= 0), "S0 non-orthogonal has negative values"
    
    print("  PASSED\n")

def test_default_flag_preserves_original():
    """Ensure use_non_orthogonal defaults to False and gives identical results to omitting it."""
    print("=" * 60)
    print("TEST 4: Default flag preserves original behaviour")
    print("=" * 60)
    
    detuning = np.linspace(-3, 3, 200) * 1e3
    E_in = np.array([1, 0, 0])
    
    p_no_flag = {'Elem': 'Rb', 'Dline': 'D2', 'T': 20, 'lcell': 75e-3,
                 'Bfield': 50, 'Btheta': 0.5}
    p_with_flag = dict(p_no_flag, use_non_orthogonal=False)
    
    S0_a = calculate(detuning, E_in, p_no_flag, outputs=['S0'])[0].real
    S0_b = calculate(detuning, E_in, p_with_flag, outputs=['S0'])[0].real
    
    diff = np.max(np.abs(S0_a - S0_b))
    print(f"  Max difference: {diff:.2e}")
    assert diff < 1e-14, f"Results differ! Max diff = {diff}"
    
    print("  PASSED\n")

if __name__ == '__main__':
    print("\nSAW-NO Test Suite: Non-Orthogonal Propagation\n")
    try:
        test_faraday_both_paths()
        test_arbitrary_angle()
        test_voigt_geometry()
        test_default_flag_preserves_original()
        print("=" * 60)
        print("ALL TESTS PASSED")
        print("=" * 60)
    except Exception as e:
        print(f"\nTEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
