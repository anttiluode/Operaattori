COMMENT
Minimal separable extracellular-field driver for Operaattori Gate 24.

Each segment receives a spatial coefficient coeff [mV]. A single GLOBAL
dimensionless waveform drive is played once. POINTER ex is connected from
Python to that segment's e_extracellular.

This follows the standard NEURON xtra pattern while keeping only the part
needed for a uniform-field assay.
ENDCOMMENT

NEURON {
    SUFFIX fielddrive
    RANGE coeff
    GLOBAL drive
    POINTER ex
}

PARAMETER {
    coeff = 0 (mV)
    drive = 0 (1)
}

ASSIGNED {
    ex (mV)
}

INITIAL {
    ex = coeff * drive
}

BEFORE BREAKPOINT {
    ex = coeff * drive
}
