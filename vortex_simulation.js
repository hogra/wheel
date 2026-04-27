/**
 * ФИЗИЧЕСКОЕ ЯДРО (Метод Роторов / Vorticity-Stream Function)
 * ИСПРАВЛЕННАЯ ВЕРСИЯ: Стабильная математика без численных взрывов
 */

const canvas = document.getElementById('vortexCanvas');
const ctx = canvas.getContext('2d');

// Параметры сетки (переведены в Grid-Units для 100% стабильности)
const res = 100;
const iter = 30; // Итерации для уравнения Пуассона
const dt = 0.2; // Внутренний шаг времени
const nu = 0.05; // Вязкость (помогает вихрям отрываться от стенки)

// Массивы сетки
let W = new Float32Array(res * res); // Завихренность (ω)
let newW = new Float32Array(res * res);
let S = new Float32Array(res * res); // Функция тока (ψ)
let U = new Float32Array(res * res); // Скорость X
let V = new Float32Array(res * res); // Скорость Y
let mask = new Float32Array(res * res); // 1 = воздух, 0 = твердое тело

// Геометрия
const cx = res / 2;
const cy = res / 2;
const radius = res * 0.45;
const innerRadius = res * 0.15;
const groundY = res * 0.85;

let currentWheelOmega = 8.0;
let numParticles = 4000;
let particles = [];

// Безопасная индексация 1D массива (исключает выход за границы)
function ix(i, j) {
    let cl_i = Math.max(0, Math.min(res - 1, Math.floor(i)));
    let cl_j = Math.max(0, Math.min(res - 1, Math.floor(j)));
    return cl_i + cl_j * res;
}

function init() {
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            let dx = i - cx;
            let dy = j - cy;
            let dist = Math.sqrt(dx * dx + dy * dy);

            if (dist < radius && dist > innerRadius && j < groundY) {
                mask[ix(i, j)] = 1.0;
            } else {
                mask[ix(i, j)] = 0.0;
            }
        }
    }

    particles = [];
    for(let k=0; k < numParticles; k++) {
        let px = Math.random() * res;
        let py = Math.random() * res;
        if (mask[ix(px, py)] === 1.0) {
            particles.push({ x: px, y: py });
        }
    }
}

// 1. Граничные условия: Генерация вихрей СТРОГО на твердой стенке
function setBoundaryVorticity() {
    for (let i = 1; i < res - 1; i++) {
        for (let j = 1; j < res - 1; j++) {
            if (mask[ix(i, j)] === 0.0) { // Мы на твердом теле
                let wall_w = 0;
                let count = 0;

                // Линейная скорость этой точки колеса
                let vx_wall = 0;
                let vy_wall = 0;
                if (j < groundY) { // Земля (j >= groundY) неподвижна
                    vx_wall = -currentWheelOmega * (j - cy) * 0.015;
                    vy_wall =  currentWheelOmega * (i - cx) * 0.015;
                }

                // ИСПРАВЛЕНИЕ: Формула Тома ДОЛЖНА иметь минус перед S (Функцией тока).
                // Это создает отрицательную обратную связь, приводя поток к равновесию.
                if (mask[ix(i+1, j)] === 1.0) { wall_w += -2 * S[ix(i+1, j)] - 2 * vy_wall; count++; }
                if (mask[ix(i-1, j)] === 1.0) { wall_w += -2 * S[ix(i-1, j)] + 2 * vy_wall; count++; }
                if (mask[ix(i, j+1)] === 1.0) { wall_w += -2 * S[ix(i, j+1)] + 2 * vx_wall; count++; }
                if (mask[ix(i, j-1)] === 1.0) { wall_w += -2 * S[ix(i, j-1)] - 2 * vx_wall; count++; }

                if (count > 0) {
                    W[ix(i, j)] = wall_w / count;
                } else {
                    W[ix(i, j)] = 0; // Глубоко внутри стенки вихрей нет
                }
            }
        }
    }
}

// 2. Уравнение Пуассона для Функции Тока
function solveStreamFunction() {
    for (let k = 0; k < iter; k++) {
        for (let i = 1; i < res - 1; i++) {
            for (let j = 1; j < res - 1; j++) {
                if (mask[ix(i, j)] === 1.0) { // Только в воздухе
                    S[ix(i, j)] = (S[ix(i+1, j)] + S[ix(i-1, j)] + S[ix(i, j+1)] + S[ix(i, j-1)] + W[ix(i, j)]) * 0.25;
                }
            }
        }
    }
}

// 3. Расчет скоростей из Функции Тока (u = dS/dy, v = -dS/dx)
function calcVelocities() {
    for (let i = 1; i < res - 1; i++) {
        for (let j = 1; j < res - 1; j++) {
            if (mask[ix(i, j)] === 1.0) {
                U[ix(i, j)] =  (S[ix(i, j+1)] - S[ix(i, j-1)]) * 0.5;
                V[ix(i, j)] = -(S[ix(i+1, j)] - S[ix(i-1, j)]) * 0.5;
            } else {
                U[ix(i, j)] = 0;
                V[ix(i, j)] = 0;
            }
        }
    }
}

// 4. Перенос вихрей воздухом (Advection) + Диффузия (Вязкость)
function advectVorticity() {
    newW.set(W);

    // Билинейная интерполяция
    function getW(x, y) {
        x = Math.max(0.5, Math.min(res - 1.5, x));
        y = Math.max(0.5, Math.min(res - 1.5, y));
        let i = Math.floor(x); let j = Math.floor(y);
        let fx = x - i; let fy = y - j;
        let w1 = W[ix(i, j)] * (1 - fx) + W[ix(i+1, j)] * fx;
        let w2 = W[ix(i, j+1)] * (1 - fx) + W[ix(i+1, j+1)] * fx;
        return w1 * (1 - fy) + w2 * fy;
    }

    for (let i = 1; i < res - 1; i++) {
        for (let j = 1; j < res - 1; j++) {
            if (mask[ix(i, j)] === 1.0) { // Двигаем вихри только внутри воздуха
                // Откуда прилетел этот воздух?
                let x = i - U[ix(i, j)] * dt;
                let y = j - V[ix(i, j)] * dt;

                // Вязкость (забирает завихренность со стенок внутрь)
                let diff = nu * dt * (W[ix(i+1, j)] + W[ix(i-1, j)] + W[ix(i, j+1)] + W[ix(i, j-1)] - 4 * W[ix(i,j)]);

                let next_w = getW(x, y) + diff;

                // ЗАЩИТА: С учетом правильной математики, жесткий clamp больше не нужен,
                // но оставим безопасный коридор от аномалий чисел с плавающей точкой.
                newW[ix(i, j)] = Math.max(-100, Math.min(100, next_w));
            }
        }
    }
    W.set(newW);
}

// 5. Движение частиц для визуализации
function advectParticles() {
    // Вспомогательная функция для плавного течения маркеров (билинейная интерполяция)
    function getInterpolatedVelocity(x, y, arr) {
        let i = Math.floor(x); let j = Math.floor(y);
        let fx = x - i; let fy = y - j;
        let v1 = arr[ix(i, j)] * (1 - fx) + arr[ix(i+1, j)] * fx;
        let v2 = arr[ix(i, j+1)] * (1 - fx) + arr[ix(i+1, j+1)] * fx;
        return v1 * (1 - fy) + v2 * fy;
    }

    for (let p of particles) {
        // ИСПРАВЛЕНИЕ: Получаем интерполированную скорость, чтобы частицы не "прыгали" по сетке
        let u = getInterpolatedVelocity(p.x, p.y, U);
        let v = getInterpolatedVelocity(p.x, p.y, V);

        p.x += u * dt;
        p.y += v * dt;

        // Жесткие рамки канваса
        if (p.x < 1) p.x = 1; if (p.x > res - 2) p.x = res - 2;
        if (p.y < 1) p.y = 1; if (p.y > res - 2) p.y = res - 2;

        // Мягкое выталкивание из стенок колеса
        if (mask[ix(p.x, p.y)] === 0.0) {
            let dx = cx - p.x; let dy = cy - p.y;
            let len = Math.sqrt(dx*dx + dy*dy) || 1;
            p.x += (dx/len) * 0.5; p.y += (dy/len) * 0.5;
        }
    }
}

function simulate() {
    setBoundaryVorticity();
    solveStreamFunction();
    calcVelocities();
    advectVorticity();
    advectParticles();
}

// --- ОТРИСОВКА ---
function getVorticityColor(w) {
    let val = Math.max(-1, Math.min(1, w * 0.4)); // 0.4 - масштаб цвета
    if (val > 0) {
        return `rgba(${Math.floor(255 * val)}, 20, 20, ${val * 0.9})`;
    } else {
        return `rgba(20, 100, ${Math.floor(255 * -val)}, ${-val * 0.9})`;
    }
}

function draw() {
    ctx.fillStyle = "#020617";
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    const cellW = canvas.width / res;
    const cellH = canvas.height / res;
    const showW = document.getElementById('showVorticity').checked;
    const showS = document.getElementById('showStreamlines').checked;
    const showP = document.getElementById('showParticles').checked;

    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            if (mask[ix(i, j)] === 0.0) {
                ctx.fillStyle = '#1e293b';
                ctx.fillRect(i * cellW, j * cellH, cellW + 1, cellH + 1);
            } else if (showW) {
                let w = W[ix(i, j)];
                if (Math.abs(w) > 0.05) {
                    ctx.fillStyle = getVorticityColor(w);
                    ctx.fillRect(i * cellW, j * cellH, cellW + 1, cellH + 1);
                }
            }
        }
    }

    if (showS) {
        ctx.fillStyle = 'rgba(255, 255, 255, 0.1)';
        for (let i = 0; i < res; i+=2) {
            for (let j = 0; j < res; j+=2) {
                if (mask[ix(i, j)] === 1.0) {
                    let s_val = Math.abs(S[ix(i, j)]);
                    if (s_val % 2.0 < 0.3) { // Линии уровня
                        ctx.fillRect(i * cellW, j * cellH, 2, 2);
                    }
                }
            }
        }
    }

    if (showP) {
        ctx.fillStyle = '#34d399';
        for (let p of particles) {
            ctx.fillRect(p.x * cellW, p.y * cellH, 1.5, 1.5);
        }
    }

    // Обводка колеса и земля
    ctx.strokeStyle = '#6ee7b7';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(cx * cellW, cy * cellH, radius * cellW, 0, 2 * Math.PI);
    ctx.stroke();

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
    currentWheelOmega = parseFloat(e.target.value);
    document.getElementById('omegaVal').innerText = currentWheelOmega.toFixed(1);
});

document.getElementById('stopBtn').addEventListener('click', () => {
    currentWheelOmega = 0.0;
    document.getElementById('omegaSlider').value = 0;
    document.getElementById('omegaVal').innerText = "0.0";
});

init();
loop();