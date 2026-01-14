let gridSize = 30;
let cells = [];
let isMouseDown = false;
let startCell = null;
let dragMode = null; // 'toggle' albo 'select'

const gridEl = document.getElementById('grid');
const outputEl = document.getElementById('output');

// Tworzenie gridu
function createGrid(n) {
    gridEl.innerHTML = '';
    gridEl.style.gridTemplateColumns = `repeat(${n}, 40px)`;
    gridEl.style.gridTemplateRows = `repeat(${n}, 40px)`;
    cells = [];

    let id = 1;
    for (let row = 0; row < n; row++) {
        cells[row] = [];
        for (let col = 0; col < n; col++) {
            const cell = document.createElement('div');
            cell.classList.add('cell');
            cell.dataset.id = id;
            cell.dataset.row = row;
            cell.dataset.col = col;
            cell.textContent = id;
            gridEl.appendChild(cell);
            cells[row][col] = cell;
            id++;
        }
    }
}

// Funkcje pomocnicze
function toggleCell(cell) {
    cell.classList.toggle('selected');
}

function selectRectangle(r1, c1, r2, c2) {
    const rowStart = Math.min(r1, r2);
    const rowEnd = Math.max(r1, r2);
    const colStart = Math.min(c1, c2);
    const colEnd = Math.max(c1, c2);
    for (let r = rowStart; r <= rowEnd; r++) {
        for (let c = colStart; c <= colEnd; c++) {
            cells[r][c].classList.add('selected');
        }
    }
}

function getSelectedIds() {
    let ids = [];
    for (let row of cells) {
        for (let cell of row) {
            if (cell.classList.contains('selected')) {
                ids.push(Number(cell.dataset.id));
            }
        }
    }
    return ids;
}

// Eventy
gridEl.addEventListener('mousedown', e => {
    if (!e.target.classList.contains('cell')) return;
    isMouseDown = true;
    startCell = e.target;
    dragMode = e.shiftKey ? 'select' : 'toggle';
    if (dragMode === 'toggle') toggleCell(e.target);
});

gridEl.addEventListener('mouseover', e => {
    if (!isMouseDown || !e.target.classList.contains('cell')) return;

    if (dragMode === 'toggle') {
        toggleCell(e.target);
    } else if (dragMode === 'select' && startCell) {
        const r1 = Number(startCell.dataset.row);
        const c1 = Number(startCell.dataset.col);
        const r2 = Number(e.target.dataset.row);
        const c2 = Number(e.target.dataset.col);

        // usuń tymczasowe zaznaczenia
        for (let row of cells) {
            for (let cell of row) {
                if (cell.classList.contains('temp-selected')) {
                    cell.classList.remove('temp-selected');
                }
            }
        }

        // zaznacz nowy prostokąt tymczasowo
        const rowStart = Math.min(r1, r2);
        const rowEnd = Math.max(r1, r2);
        const colStart = Math.min(c1, c2);
        const colEnd = Math.max(c1, c2);
        for (let r = rowStart; r <= rowEnd; r++) {
            for (let c = colStart; c <= colEnd; c++) {
                cells[r][c].classList.add('temp-selected');
            }
        }
    }
});

document.addEventListener('mouseup', () => {
    if (dragMode === 'select') {
        for (let row of cells) {
            for (let cell of row) {
                if (cell.classList.contains('temp-selected')) {
                    cell.classList.remove('temp-selected');
                    cell.classList.add('selected');
                }
            }
        }
    }
    isMouseDown = false;
    startCell = null;
    dragMode = null;
});

// Shift+klik prostokąt
gridEl.addEventListener('click', e => {
    if (!e.target.classList.contains('cell')) return;
    if (e.shiftKey) {
        if (!startCell) {
            startCell = e.target;
        } else {
            const r1 = Number(startCell.dataset.row);
            const c1 = Number(startCell.dataset.col);
            const r2 = Number(e.target.dataset.row);
            const c2 = Number(e.target.dataset.col);
            selectRectangle(r1, c1, r2, c2);
            startCell = null;
        }
    }
});

// Generowanie tablicy ID
document.getElementById('generateArray').addEventListener('click', () => {
    const ids = getSelectedIds();
    outputEl.textContent = JSON.stringify(ids);
});

// Generowanie nowego gridu
document.getElementById('generateGrid').addEventListener('click', () => {
    gridSize = Number(document.getElementById('size').value);
    createGrid(gridSize);
});

// Tworzymy domyślny grid
createGrid(gridSize);
