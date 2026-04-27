const canvas = document.getElementById('simCanvas');
const ctx = canvas.getContext('2d');

const res = 100;
const iter = 40;
const simHeight = 1.0;
const h = simHeight / res;
const dt = 1.0 / 60.0;

let U = new Float32Array((res + 1) * res);
let V = new Float32Array(res * (res + 1));
let newU = new Float32Array((res + 1) * res);
let newV = new Float32Array(res * (res + 1));
let P = new Float32Array(res * res);
let s = new Float32Array(res * res);

const cx = res / 2;
const cy = res / 2;
const radius = res * 0.45;
const innerRadius = res * 0.15;
const groundY = res * 0.9;

let currentOmega = 5.0;
let numParticles = 8000;
let particles = [];

function ixU(i, j) { return i + j * (res + 1); }
function ixV(i, j) { return i + j * res; }
function ix(i, j) { return i + j * res; }

function init() {
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            let dx = i - cx;
            let dy = j - cy;
            let dist = Math.sqrt(dx * dx + dy * dy);

            if (dist < radius && dist > innerRadius && j < groundY) {
                s[ix(i, j)] = 1.0;
            } else {
                s[ix(i, j)] = 0.0;
            }
        }
    }

    particles = [];
    while (particles.length < numParticles) {
        let px = Math.random() * res;
        let py = Math.random() * res;
        let i = Math.floor(px);
        let j = Math.floor(py);
        if (i >= 0 && i < res && j >= 0 && j < res && s[ix(i, j)] === 1.0) {
            particles.push({ x: px, y: py });
        }
    }
}

function setBoundaryVelocities() {
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            if (s[ix(i, j)] === 0.0) {
                let vx = 0, vy = 0;
                if (j < groundY) {
                    vx = -currentOmega * (j - cy) * h;
                    vy =  currentOmega * (i - cx) * h;
                }
                if (i < res - 1 && s[ix(i + 1, j)] === 1.0) U[ixU(i + 1, j)] = vx;
                if (i > 0 && s[ix(i - 1, j)] === 1.0) U[ixU(i, j)] = vx;
                if (j < res - 1 && s[ix(i, j + 1)] === 1.0) V[ixV(i, j + 1)] = vy;
                if (j > 0 && s[ix(i, j - 1)] === 1.0) V[ixV(i, j)] = vy;
            }
        }
    }
}

function solveIncompressibility() {
    for (let k = 0; k < iter; k++) {
        for (let i = 1; i < res - 1; i++) {
            for (let j = 1; j < res - 1; j++) {
                if (s[ix(i, j)] === 0.0) continue;

                let sx0 = s[ix(i - 1, j)];
                let sx1 = s[ix(i + 1, j)];
                let sy0 = s[ix(i, j - 1)];
                let sy1 = s[ix(i, j + 1)];
                let s_sum = sx0 + sx1 + sy0 + sy1;

                if (s_sum === 0.0) continue;

                let div = U[ixU(i + 1, j)] - U[ixU(i, j)] + V[ixV(i, j + 1)] - V[ixV(i, j)];
                let p_val = -div / s_sum;
                p_val *= 1.9;

                P[ix(i, j)] += p_val;

                if (sx0) U[ixU(i, j)] -= p_val;
                if (sx1) U[ixU(i + 1, j)] += p_val;
                if (sy0) V[ixV(i, j)] -= p_val;
                if (sy1) V[ixV(i, j + 1)] += p_val;
            }
        }
    }
}

function advectVelocities() {
    newU.set(U);
    newV.set(V);

    const getU = (x, y) => {
        x = Math.max(0, Math.min(x, res));
        y = Math.max(0.5, Math.min(y, res - 0.5));
        let i = Math.floor(x), j = Math.floor(y - 0.5);
        let fx = x - i, fy = y - 0.5 - j;
        let u1 = U[ixU(i, j)] * (1 - fx) + U[ixU(Math.min(i + 1, res), j)] * fx;
        let u2 = U[ixU(i, Math.min(j + 1, res - 1))] * (1 - fx) + U[ixU(Math.min(i + 1, res), Math.min(j + 1, res - 1))] * fx;
        return u1 * (1 - fy) + u2 * fy;
    };

    const getV = (x, y) => {
        x = Math.max(0.5, Math.min(x, res - 0.5));
        y = Math.max(0, Math.min(y, res));
        let i = Math.floor(x - 0.5), j = Math.floor(y);
        let fx = x - 0.5 - i, fy = y - j;
        let v1 = V[ixV(i, j)] * (1 - fy) + V[ixV(i, Math.min(j + 1, res))] * fy;
        let v2 = V[ixV(Math.min(i + 1, res - 1), j)] * (1 - fy) + V[ixV(Math.min(i + 1, res - 1), Math.min(j + 1, res))] * fy;
        return v1 * (1 - fx) + v2 * fx;
    };

    for (let i = 1; i < res; i++) {
        for (let j = 1; j < res - 1; j++) {
            if (s[ix(i, j)] !== 0.0 && s[ix(i - 1, j)] !== 0.0) {
                let x = i, y = j + 0.5;
                let u = U[ixU(i, j)];
                let v = (V[ixV(i - 1, j)] + V[ixV(i, j)] + V[ixV(i - 1, j + 1)] + V[ixV(i, j + 1)]) * 0.25;
                x -= u * dt / h; y -= v * dt / h;
                newU[ixU(i, j)] = getU(x, y);
            }
        }
    }

    for (let i = 1; i < res - 1; i++) {
        for (let j = 1; j < res; j++) {
            if (s[ix(i, j)] !== 0.0 && s[ix(i, j - 1)] !== 0.0) {
                let x = i + 0.5, y = j;
                let u = (U[ixU(i, j - 1)] + U[ixU(i, j)] + U[ixU(i + 1, j - 1)] + U[ixU(i + 1, j)]) * 0.25;
                let v = V[ixV(i, j)];
                x -= u * dt / h; y -= v * dt / h;
                newV[ixV(i, j)] = getV(x, y);
            }
        }
    }
    U.set(newU);
    V.set(newV);
}

function advectParticles() {
    for (let p of particles) {
        let i = Math.max(0, Math.min(Math.floor(p.x), res - 1));
        let j = Math.max(0, Math.min(Math.floor(p.y), res - 1));
        let u = (U[ixU(i, j)] + U[ixU(i + 1, j)]) * 0.5;
        let v = (V[ixV(i, j)] + V[ixV(i, j + 1)]) * 0.5;
        p.x += u * dt / h;
        p.y += v * dt / h;

        if (p.x < 0) p.x = 0; if (p.x > res) p.x = res;
        if (p.y < 0) p.y = 0; if (p.y > res) p.y = res;

        if (s[ix(Math.floor(p.x), Math.floor(p.y))] === 0.0) {
            let dx = cx - p.x, dy = cy - p.y;
            let len = Math.sqrt(dx*dx + dy*dy) || 1;
            p.x += (dx/len) * 0.5; p.y += (dy/len) * 0.5;
        }
    }
}

function simulate() {
    setBoundaryVelocities();
    advectVelocities();
    setBoundaryVelocities();
    solveIncompressibility();
    setBoundaryVelocities();
    advectParticles();
}

function getSciColor(val, minVal, maxVal) {
    let t = Math.max(0.0, Math.min(1.0, (val - minVal) / (maxVal - minVal)));
    let r = Math.floor(Math.max(0, 255 * (2 * t - 1)));
    let g = Math.floor(Math.max(0, 255 * (1 - Math.abs(2 * t - 1))));
    let b = Math.floor(Math.max(0, 255 * (1 - 2 * t)));
    // для прозрачности фона давления
    return `rgba(${r},${g},${b}, 0.5)`;
}

function draw() {
    // Темный фон
    ctx.fillStyle = "#020617";
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    const cellW = canvas.width / res;
    const cellH = canvas.height / res;
    const showP = document.getElementById('showPressure').checked;
    const showV = document.getElementById('showVelocity').checked;
    const showM = document.getElementById('showParticles').checked;

    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            if (s[ix(i, j)] === 0.0) {
                ctx.fillStyle = '#1e293b';
                ctx.fillRect(i * cellW, j * cellH, cellW + 1, cellH + 1);
            } else if (showP) {
                ctx.fillStyle = getSciColor(P[ix(i, j)], -0.05, 0.05);
                ctx.fillRect(i * cellW, j * cellH, cellW + 1, cellH + 1);
            }
        }
    }

    if (showM) {
        ctx.fillStyle = '#22d3ee';
        for (let p of particles) {
            ctx.beginPath();
            // для скорости
            ctx.fillRect(p.x * cellW - 0.5, p.y * cellH - 0.5, 1.5, 1.5);
        }
    }

    if (showV) {
        ctx.strokeStyle = 'rgba(253, 224, 71, 0.5)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let i = 0; i < res; i += 4) {
            for (let j = 0; j < res; j += 4) {
                if (s[ix(i, j)] === 1.0) {
                    let u = (U[ixU(i, j)] + U[ixU(i + 1, j)]) * 0.5;
                    let v = (V[ixV(i, j)] + V[ixV(i, j + 1)]) * 0.5;
                    ctx.moveTo((i + 0.5) * cellW, (j + 0.5) * cellH);
                    ctx.lineTo((i + 0.5) * cellW + u * 150, (j + 0.5) * cellH + v * 150);
                }
            }
        }
        ctx.stroke();
    }

    ctx.strokeStyle = '#3b82f6';
    ctx.setLineDash([5, 5]); // для эстетики
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.arc(cx * cellW, cy * cellH, radius * cellW, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.strokeStyle = '#ef4444';
    ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.moveTo(0, groundY * cellH);
    ctx.lineTo(canvas.width, groundY * cellH);
    ctx.stroke();
}

function loop() {
    simulate();
    draw();
    requestAnimationFrame(loop);
}

document.getElementById('omegaSlider').addEventListener('input', (e) => {
    currentOmega = parseFloat(e.target.value);
    document.getElementById('omegaVal').innerText = currentOmega.toFixed(1);
});

document.getElementById('stopBtn').addEventListener('click', () => {
    currentOmega = 0.0;
    document.getElementById('omegaSlider').value = 0;
    document.getElementById('omegaVal').innerText = "0.0";
});

init();
loop();