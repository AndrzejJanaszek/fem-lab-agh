# AGH Labolatoria Metoda Elementów Skończonych (FEM)

[> English version](./README.md) 

Program MES do symulacji procesów cieplnych w 2 wymiarach. 
Bazowa wersja programu zakłada możliwość symulacji stacjonarnego i niestacjonarnego procesu dla jednego typu materiału z zadanymi warunkami brzegowymi.
Wersja dla problemu rzeczywistego rozbudowuje program o możliwość dodania drugiego materiału, dodania źródła ciepła oraz możliwość zastosowania sztucznej dyfuzji*.

*Sztuczna dyfuzja* - uśrednianie temperatury w każdym kroku symulacji w zadanej obiętości. W rozpatrywanym problemie rzeczywistym jest to uśrednianie temperatury w całej obiętości powietrza w siatce MES.

Program zapisuje dane symulacji do plików `.pvd`/`.vtu`, aby potem można było przeanalizować symulację w czasie i ją zwizualizować.

Dla symulacji nieustalonej program zatrzymuje się gdy zostanie usiągnięta temperatura maksymalna lub proces stanie się ustalony (zmiana maksymalnej temperatury w miedzy krokami symulacji będzie poniżej zadanego epsilona).

### Rozkłady temperatur dla procesu ustalonego w siatkach testowych (wersja bazowa programu) wraz z schematem siatek.

<div style="display: flex; justify-content: center; gap: 10px;">
  <img src="./others/sprawozdanie/img/4_4_siatka_schemat.png" width="30%">
  <img src="./others/sprawozdanie/img/4_4_mix_siatka_schemat.png" width="30%">
  <img src="./others/sprawozdanie/img/30x30_siatka_schemat.png" width="30%">
</div>

### Format pliku danych siatki (z definicjami i przykładami dla bardziej złożonych wartości)
```
SimulationTime <czas symulacji w sekundach>
SimulationStepTime <krok czasowy symulacji w sekundach>
Conductivity <przewodność materiału>
Alfa <współczynnik konwekcyjnej wymiany ciepła>
Tot <temperatura otoczenia>
InitialTemp <początkowa temperatura materiału>
Density <gęstość materiału>
SpecificHeat <ciepło właściwe materiału>
Nodes_number <liczba węzłów w siatce MES>
Elements_number <liczba elementów w siatce>
*Node
<id węzła (1...)>, <współrzędna x węzła>, <współrzędna y węzła>
      1,  0.100000001, 0.00499999989
                    ...
     16,           0., -0.0949999988
*Element, type=DC2D4
<id elementu (1...)>, <id węzła>, <id węzła>, <id węzła>, <id węzła>
 1,  1,  2,  6,  5
        ...
 9, 11, 12, 16, 15
*BC
<id węzłów oddzielone przecinkiem (,) – krawędzie z warunkiem brzegowym>
1, 2, ... 15, 16
```

### Właściwości fizyczne materiałów i konfiguracja dla problemu rzeczywistego

Z powodu szybkiego i dość chaotycznego podejścia do rozwiązania problemu rzeczywistego, właściwości materiałowe są zapisane na stałe w kodzie, jak pokazano poniżej (wersja uproszczona).

``` cpp
const double COPPER_CONDUCTIVITY = 394.85;
const double COPPER_DENSITY = 8911.47;
const double COPPER_SPECIFIC_HEAT = 384.37;
...
const std::vector<int> COPPER_ELEMENTS_2 = {95,96 ... 900};
const std::vector<int> COPPER_ELEMENTS_1 = {100,101, ... 900};

const std::vector<int> INITIAL_HOT_ELEMENTS = {881,882 ... 890};
```

### Opis problemu rzeczywistego

Problem rzeczywisty zakłada przeprowadzenie symulacji dla radiatora procesora. Porównanie temperatur (max, min, avg) w zależności od liczby żeber, oraz porównanie symulacji dla procesu ustalonego, nieustalonego oraz nieustalonego ze sztuczną dyfucją.

<img src="./others/sprawozdanie/img/schemat_radiatora.png" width="400">

Symulacjca przeprowadzana jest dla siatki 30x30 (elementów). Mimo teoretycznego rozróżnienia na *żebra*, *IHS* oraz *cpu die*, wszystkie te elementy są po prostu meteriałem z miedzi.
Elementy siatki *cpu die* mają zadane źródło ciepła.

Siatka jest przeskalowana w taki spsób aby wymiary *cpu die* oraz *IHS* były jak najbliższe rzeczywistym.
*IHS* ma zazwyczaj wymiary 40x40mm, a *cpu die* 13x13.
W symulacji długość boku elementu to 1.3mm, dzięki czemu długość *cpu die* to 13mm a *IHS* 39mm.

#### Schematy siatki dla 3 scenariuszy (1, 2 i 4 żeber radiatora).
<div style="display: flex; justify-content: center; gap: 10px;">
  <img src="./others/sprawozdanie/img/schemat_1k_radiator.png" width="30%">
  <img src="./others/sprawozdanie/img/schemat_2k_radiator.png" width="30%">
  <img src="./others/sprawozdanie/img/schemat_4k_radiator.png" width="30%">
</div>

#### Rozkład temperatur dla procesu ustalonego
<div style="display: flex; justify-content: center; gap: 10px;">
  <img src="./others/sprawozdanie/img/radiator_1zebra_ustalony.png" width="30%">
  <img src="./others/sprawozdanie/img/radiator_2zebra_ustalony.png" width="30%">
  <img src="./others/sprawozdanie/img/radiator_4zebra_ustalony.png" width="30%">
</div>

---

### Dodatkowe materiały

Więcej informacji można znaleźć w sprawozdaniu (*wersja polska*):  
[Sprawozdanie (PDF)](./others/sprawozdanie/sprawozdanie_mes_andrzej_janaszek_16_01_2026.pdf)  

Lub bezpośrednio pod ścieżką: `./others/sprawozdanie/sprawozdanie_mes_andrzej_janaszek_16_01_2026.pdf`