fn efield(t) {
    let hbar  = 0.658212;           // reduced planck constant (eV/fs)
    let amp   = 0.005;              // amplitude = 0.005 (Volt/Angstrom)
    let hnu   = 1.71;             // photon energy = h * nu = 1.4 (eV)
    let omega = hnu / hbar;         // omega = hnu / hbar (rad/fs)
    let duration = 1000;
    let omega_envelop = 2.0 / duration * pi;
    
    let x = amp * cos(omega * t);   // electric field at `x` direction at time `t`
    let y = amp * sin(omega * t);   // electric field at `y` direction at time `t`
    let z = 0.0;                    // no electric field at `z` direction

    let envelop = if t <= duration {
        sin(omega_envelop*t - 0.5*pi) * 0.5 + 0.5;
    } else {
        0
    };

    return [x*envelop, y*envelop, z*envelop];               // this statement is required.
}