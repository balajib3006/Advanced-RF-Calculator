// Tab switching
document.getElementById('tabPatchBtn').addEventListener('click', ()=>switchTab('patch'));
document.getElementById('tabLinkBtn').addEventListener('click', ()=>switchTab('link'));
function switchTab(name){
    document.getElementById('tab-patch').style.display = name==='patch' ? '' : 'none';
    document.getElementById('tab-link').style.display  = name==='link' ? '' : 'none';
    document.querySelectorAll('.tab-btn').forEach(b=>b.classList.remove('active'));
    (name==='patch' ? document.getElementById('tabPatchBtn') : document.getElementById('tabLinkBtn')).classList.add('active');
}

  // Helpers
  const c = 299792458;
  const fmt = (v,d=3)=>Number.isFinite(v)?Number(v).toFixed(d):'—';
  const fmt2 = (v,d=2)=>Number.isFinite(v)?Number(v).toFixed(d):'—';
  const clamp = (x,lo,hi)=>Math.min(Math.max(x,lo),hi);

  // Modulation presets
  function applyModPreset(){
    const sel = document.getElementById('l_mod').value;
    if(sel!==''){ document.getElementById('l_reqSNR').value = parseFloat(sel); }
  }

  // Microstrip width (Hammerstad/Jensen)
  function microstripWidth(Z0, er, h){
    const A = (Z0/60)*Math.sqrt((er+1)/2) + ((er-1)/(er+1))*(0.23 + 0.11/er);
    let w_h;
    if (Z0 <= (120*Math.PI)/Math.sqrt(er+1)){
      w_h = (8*Math.exp(A)) / (Math.exp(2*A) - 2);
    } else {
      const B = (377*Math.PI) / (2*Z0*Math.sqrt(er));
      w_h = (2/Math.PI) * ( B - 1 - Math.log(2*B - 1) + ((er-1)/(2*er)) * (Math.log(B-1) + 0.39 - 0.61/er) );
    }
    const Wf = w_h * h;
    const eps_eff_line = (er+1)/2 + (er-1)/2 * (1/Math.sqrt(1 + 12*h/Wf));
    return {Wf, eps_eff_line, w_h};
  }

  /* ---------- Patch calculation ---------- */
  function clearPatch() {
 const f = document.getElementById('patchForm');
 if (f) f.reset();
 document.getElementById('patchResults').innerHTML = '';
}
  function calcPatch(){
    const fGHz = parseFloat(document.getElementById('p_freq').value);
    const er = parseFloat(document.getElementById('p_er').value);
    const h_mm = parseFloat(document.getElementById('p_hmm').value);
    const tanD = parseFloat(document.getElementById('p_tanD').value);
    const Z0 = parseFloat(document.getElementById('p_Z0').value);
    const eta_ap = parseFloat(document.getElementById('p_etaap').value);
    if (![fGHz,er,h_mm,tanD,Z0,eta_ap].every(v=>Number.isFinite(v))){ alert('Fill all patch fields'); return; }
    if (!(fGHz>0 && er>1 && h_mm>0 && tanD>0 && Z0>0 && eta_ap>0 && eta_ap<=1)){ alert('Check inputs'); return; }
    const f = fGHz*1e9; const h = h_mm/1000; const lambda0 = c/f;
    const W = (c/(2*f))*Math.sqrt(2/(er+1));
    const erEff = (er+1)/2 + (er-1)/2 * (1/Math.sqrt(1 + 12*h/W));
    const W_h = W/h;
    const dL = 0.412 * h * ((erEff + 0.3)*(W_h + 0.264)) / ((erEff - 0.258)*(W_h + 0.8));
    const Leff = c / (2 * f * Math.sqrt(erEff));
    const L = Leff - 2*dL; if(!(L>0)){ alert('Non-physical L'); return; }
    const Lsub = 6*h+L, Wsub = 6*h+W;
    const lambdaEff = lambda0 / Math.sqrt(erEff);
    const Qd = 1/tanD; const Qr = (Math.PI*Math.sqrt(erEff)*L)/(2*h); const Q_total = (Qd*Qr)/(Qd+Qr);
    const BWpct = 3.771 * (er - 1 / Math.pow(er, 2)) * (h/lambda0) * (W / L);
    const Rin_edge = 90 * (Math.pow(er, 2) / er - 1) * (L / W);
    let yin_m=0, inset_note='';
    if(Z0 <= Rin_edge){ yin_m = (L / Math.PI) * Math.acos(Math.sqrt(clamp(Z0/Rin_edge,0,1))); } else { yin_m=0; inset_note=' (Z₀ > R_edge — inset matching not possible)'; }
    const insetLen = W/2, insetW = W/10;
    function fmn(m,n){ return (c/(2*Math.sqrt(erEff))) * Math.sqrt(Math.pow(m/L,2)+Math.pow(n/W,2)); }
    const modes = [{name:'TM10 (fund.)',f:fmn(1,0)},{name:'TM01',f:fmn(0,1)},{name:'TM20',f:fmn(2,0)},{name:'TM11',f:fmn(1,1)},{name:'TM02',f:fmn(0,2)}];
    const Aphys = L*W; const D_lin = eta_ap * (4*Math.PI*Aphys)/(lambda0*lambda0); const D_dBi = 10*Math.log10(D_lin); const eta_rad = Qr/(Qr+Qd); const G_lin = eta_rad * D_lin; const G_dBi = 10*Math.log10(G_lin);
    const feed = microstripWidth(Z0, er, h); const Wf_mm = feed.Wf*1000;
    // render
    const container = document.getElementById('patchResults'); container.innerHTML='';
    function makeCard(title, rowsHtml){ const d=document.createElement('div'); d.className='card'; d.innerHTML=`<h3>${title}</h3>${rowsHtml}`; return d; }
    function r(n,v){ return `<div class="row"><div class="name">${n}</div><div class="val">${v}</div></div>`; }
    let html='';
    html += r('Wavelength (λ₀)', `${fmt2(lambda0*1000,2)} mm`);
    html += r('Effective Wavelength', `${fmt2(lambdaEff*1000,2)} mm`);
    html += r('Effective Permittivity ε_eff', `${fmt(erEff,4)}`);
    container.appendChild(makeCard('Wavelengths & Permittivity', html));
    html=''; html+=r('Patch Width W', `${fmt2(W*1000,3)} mm`); html+=r('Patch Length L', `${fmt2(L*1000,3)} mm`); html+=r('Effective Length L_eff', `${fmt2(Leff*1000,3)} mm`); container.appendChild(makeCard('Patch Geometry', html));
    html=''; html+=r('Substrate Length (≈6×)', `${fmt2(Lsub*1000,3)} mm`); html+=r('Substrate Width (≈6×)', `${fmt2(Wsub*1000,3)} mm`); html+=`<div class="note">Use ≥6× patch dims to reduce edge effects; increase for thick/low-εr substrates.</div>`; container.appendChild(makeCard('Substrate (rule-of-thumb)', html));
    html=''; html+=r('Dielectric Q (Qd)', `${fmt(Qd,2)}`); html+=r('Radiation Q (Qr)', `${fmt(Qr,2)}`); html+=r('Combined Q', `${fmt(Q_total,2)}`); html+=r('Bandwidth (est.)', `${fmt(BWpct,3)}`); html+=`<div class="note">BW & Q do not include conductor/surface-wave losses. Use EM simulation for final BW.</div>`; container.appendChild(makeCard('Q & Bandwidth', html));
    html=''; html+=r('Edge Input Resistance R_edge (est.)', `${fmt(Rin_edge,2)} Ω`); html+=r(`Inset Depth y (Z₀=${fmt2(Z0,0)} Ω)`, `${fmt2(yin_m*1000,3)} mm${inset_note}`); html+=r('Inset Feed Length (suggest)', `${fmt2(insetLen*1000,3)} mm`); html+=r('Inset Feed Width (suggest)', `${fmt2(insetW*1000,3)} mm`); html+=`<div class="note">Model: R_in(y)=R_edge ⋅ cos²(π y / L). If Z₀ > R_edge, inset matching is not possible.</div>`; container.appendChild(makeCard('Inset Feed & Edge Resistance', html));
    html=''; html+=r(`Microstrip Feedline W_f (Z₀=${fmt2(Z0,0)} Ω)`, `${fmt2(Wf_mm,3)} mm`); html+=r('Feedline w/h', `${fmt(feed.w_h,4)}`); html+=r('Feedline ε_eff', `${fmt(feed.eps_eff_line,4)}`); html+=`<div class="note">Width computed with Hammerstad/Jensen closed-form.</div>`; container.appendChild(makeCard('Microstrip Feedline', html));
    html=''; modes.forEach(m=>html+=r(m.name,`${fmt2(m.f/1e9,4)} GHz`)); html+=`<div class="note">Cavity-model estimates; verify with full-wave EM.</div>`; container.appendChild(makeCard('Multi-mode Frequencies', html));
    html=''; html+=r('Aperture Area A', `${fmt2(Aphys*1e6,2)} mm²`); html+=r('Aperture Efficiency η_ap (user)', `${fmt(eta_ap,3)}`); html+=r('Radiation Efficiency η_rad (≈Qr/(Qr+Qd))', `${fmt(eta_rad,3)}`); html+=r('Directivity (dBi)', `${fmt2(D_dBi,2)} dBi`); html+=r('Gain (dBi, est.)', `${fmt2(G_dBi,2)} dBi`); html+=`<div class="note">Gain uses η_rad from Qs; conductor & surface waves not included.</div>`; container.appendChild(makeCard('Radiation Performance (est.)', html));
  }

  /* ---------- Link budget + SNR + sensitivity + NF calculations ---------- */
  function clearLink() {
 const f = document.getElementById('linkForm');
 if (f) f.reset();
 document.getElementById('linkResults').innerHTML = '';
}

function calcLink() {
  // === Grab Inputs ===
  const freqGHz = parseFloat(document.getElementById("l_freq").value) || 2.4; // GHz
  const dist_m  = parseFloat(document.getElementById("l_dist").value) || 1000;

  // Tx side
  let txVal = parseFloat(document.getElementById("l_txval").value) || 0;
  const txUnit = document.getElementById("l_txunit").value;
  if (txUnit === "dBW") txVal = txVal + 30; // convert to dBm
  const Gt = parseFloat(document.getElementById("l_Gt").value) || 0; // dBi
  const txLoss = parseFloat(document.getElementById("l_txLoss")?.value) || 0;
  const Nt = parseInt(document.getElementById("l_Nt")?.value) || 1;
  const eta = parseFloat(document.getElementById("l_eta")?.value) || 1;
  const polLoss = parseFloat(document.getElementById("l_polLoss")?.value) || 0;
  const pointingLoss = parseFloat(document.getElementById("l_pointing")?.value) || 0;
  const backoff = parseFloat(document.getElementById("l_backoff")?.value) || 0;
  const EIRP_limit = parseFloat(document.getElementById("l_EIRPmax")?.value) || null;

  // Channel / Propagation
  const model = document.getElementById("l_model")?.value || "FSPL";
  const atmLoss = parseFloat(document.getElementById("l_atmLoss")?.value) || 0;
  const rainLoss = parseFloat(document.getElementById("l_rainLoss")?.value) || 0;
  const folLoss  = parseFloat(document.getElementById("l_folLoss")?.value) || 0;
  const clutter  = parseFloat(document.getElementById("l_clutter")?.value) || 0;

  // Rx side
  const Gr = parseFloat(document.getElementById("l_Gr").value) || 0;
  const rxLoss = parseFloat(document.getElementById("l_rxLoss")?.value) || 0;
  const Nr = parseInt(document.getElementById("l_Nr")?.value) || 1;
  const NF = parseFloat(document.getElementById("l_NF").value) || 5;

  // System / Modulation
  const BW = parseFloat(document.getElementById("l_BW").value) || 1e6; // Hz
  const reqSNR = parseFloat(document.getElementById("l_reqSNR").value) || 10;

  // === Normalization ===
  const freqHz = freqGHz * 1e9;
  const dist_km = dist_m / 1000;
  const realizedGainTx = Gt + 10 * Math.log10(eta);

  // === Step 1: EIRP ===
  let EIRP = txVal - txLoss + realizedGainTx - polLoss - pointingLoss - backoff;
  // Add MIMO Tx array gain if beamforming
  if (Nt > 1) {
    EIRP += 10 * Math.log10(Nt); // ideal BF
  }
  if (EIRP_limit && EIRP > EIRP_limit) {
    EIRP = EIRP_limit; // regulatory cap
  }

  // === Step 2: Path Loss (select model) ===
  let pathLoss = 0;
  if (model === "FSPL") {
    pathLoss = 32.45 + 20 * Math.log10(freqGHz*1000) + 20 * Math.log10(dist_km); // freq in MHz
  } else if (model === "LogDistance") {
    const n = parseFloat(document.getElementById("l_pathExp")?.value) || 3;
    const d0 = parseFloat(document.getElementById("l_d0")?.value) || 1;
    const L0 = 32.45 + 20 * Math.log10(freqGHz*1000) + 20 * Math.log10(d0/1000);
    pathLoss = L0 + 10 * n * Math.log10(dist_m/d0);
  } else {
    // TODO: implement Hata, COST231, 3GPP, ITU, Two-ray etc.
    pathLoss = 32.45 + 20 * Math.log10(freqGHz*1000) + 20 * Math.log10(dist_km);
  }

  // === Step 3: Total Channel Loss ===
  const channelLoss = pathLoss + atmLoss + rainLoss + folLoss + clutter;

  // === Step 4: Rx Power ===
  let Pr = EIRP - channelLoss - rxLoss + Gr;
  if (Nr > 1) {
    Pr += 10 * Math.log10(Nr); // Rx array gain
  }

  // === Step 5: Noise Power ===
  const N0 = -174 + 10 * Math.log10(BW) + NF; // dBm

  // === Step 6: SNR ===
  let SNR = Pr - N0;

  // === Step 7: Sensitivity & Margin ===
  const sens = N0 + reqSNR;
  const linkMargin = Pr - sens;

  // === Step 8: Capacity ===
  const Cbps = BW * Math.log2(1 + Math.pow(10, SNR/10));
  const Cmbps = Cbps / 1e6;

  // === Results ===
  let html = `
    <div class="card"><h3>Transmitter</h3>
      <p>EIRP: ${EIRP.toFixed(2)} dBm</p>
      <p>Tx Antenna Realized Gain: ${realizedGainTx.toFixed(2)} dBi</p>
    </div>
    <div class="card"><h3>Channel</h3>
      <p>Path Loss (${model}): ${pathLoss.toFixed(2)} dB</p>
      <p>Total Channel Loss: ${channelLoss.toFixed(2)} dB</p>
    </div>
    <div class="card"><h3>Receiver</h3>
      <p>Received Power: ${Pr.toFixed(2)} dBm</p>
      <p>Noise Power: ${N0.toFixed(2)} dBm</p>
      <p>SNR: ${SNR.toFixed(2)} dB</p>
      <p>Receiver Sensitivity: ${sens.toFixed(2)} dBm</p>
      <p>Link Margin: ${linkMargin.toFixed(2)} dB</p>
    </div>
    <div class="card"><h3>System Capacity</h3>
      <p>Shannon Capacity: ${Cmbps.toFixed(2)} Mbps</p>
    </div>
  `;

  document.getElementById("linkResults").innerHTML = html;
}

  // start on patch
  switchTab('patch');