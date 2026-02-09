let GRID_NODES = 31;
let GRID_ELEMENTS = 30;
const grid = document.getElementById("grid");
const nodes = [];

function buildGrid(zeroBased, showElements) {
    grid.innerHTML = "";
    nodes.length = 0;

    const offset = zeroBased ? 0 : 1;

    const size = showElements ? GRID_ELEMENTS : GRID_NODES;

    // ustawienie siatki CSS
    grid.style.gridTemplateColumns = `repeat(${size}, 16px)`;
    grid.style.gridTemplateRows = `repeat(${size}, 16px)`;

    for (let row = 0; row < size; row++) {
        for (let col = 0; col < size; col++) {
            const id = row * size + (size - col - 1) + offset; // od prawego górnego rogu
            const div = document.createElement("div");
            div.classList.add("node");
            div.dataset.id = id;
            grid.appendChild(div);
            nodes.push(div);
        }
    }
}

// initial build (nody)
buildGrid(false, false);

function highlightNodes() {
    const zeroBased = document.getElementById("zeroBased").checked;
    const showElements = document.getElementById("showElements").checked;

    buildGrid(zeroBased, showElements);

    nodes.forEach(n => n.classList.remove("active"));

    const singleId = document.getElementById("nodeId").value;
    const listInput = document.getElementById("nodeList").value;

    const ids = new Set();

    if (singleId !== "") ids.add(parseInt(singleId));
    if (listInput !== "") {
        listInput.split(",")
            .map(v => v.trim())
            .filter(v => v !== "")
            .forEach(v => ids.add(parseInt(v)));
    }

    nodes.forEach(n => {
        if (ids.has(parseInt(n.dataset.id))) n.classList.add("active");
    });
}
