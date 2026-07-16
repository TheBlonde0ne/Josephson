# Auswertung: Magnetfeldabhängigkeit des Josephson-Stroms

Fortführung der Auswertung und physikalische Diskussion

Die Anpassung des theoretischen Fraunhofer-Musters gemäß der Gleichung

$$I_c(B) = I_0 \left\vert{} \text{sinc}\left( \frac{\pi (B - B_{\text{offset}})}{\Delta B} \right) \right\vert{} + C$$

an die experimentellen Daten liefert die in Tabelle \ref{tab:josephson_fit_ergebnisse} zusammengefassten Parameter. Im Folgenden werden diese physikalisch interpretiert, mit den mikroskopischen suparleitenden Längenskalen verknüpft und im Kontext der experimentellen Rahmenbedingungen diskutiert.

1. Maximaler Josephson-Strom $I_0$

Der aus dem Fit extrahierte maximale Josephson-Strom bei verschwindendem effektiven Magnetfeld beträgt $I_0 = (59{,}09 \pm 3{,}53)\ \mu\text{A}$ (bzw. $0{,}0591 \pm 0{,}0035\ \text{mA}$). Dieser Parameter beschreibt die maximale Paarstrom-Kopplungsstärke des schwach gekoppelten Supraleiters (Weak Link) im feldfreien Zustand. Er ist direkt proportional zur temperaturabhängigen Energielücke $\Delta(T)$ des Supraleiters und umgekehrt proportional zum Normalleitungswiderstand $R_N$ des Übergangs (Ambegaokar-Baratoff-Relation). Die statistische Unsicherheit von ca. $6\ \%$ deutet auf eine robuste Datengrundlage im Bereich des Hauptmaximums hin.

2. Magnetisches Feld-Offset $B_{\text{offset}}$

Für das magnetische Feld-Offset ergibt sich ein Wert von $B_{\text{offset}} = (-5{,}2 \pm 13{,}7)\ \mu\text{T}$ (bzw. $-0{,}0052 \pm 0{,}0137\ \text{mT}$). Ein solches Feld-Offset resultiert physikalisch aus zwei primären Quellen:

Dem ungeschirmten Anteil des externen Erdmagnetfeldes (typischerweise im Bereich von $30\ \mu\text{T}$ bis $60\ \mu\text{T}$),

Remanenten Magnetisierungen der im Versuchsaufbau verwendeten Materialien nahe der Probe.

Bemerkenswert ist hierbei, dass die statistische Unsicherheit ($\sigma = 13{,}7\ \mu\text{T}$) die absolute Magnitude des Offsets übersteigt. Das extrahierte Feld-Offset ist somit im Rahmen der $1\sigma$-Standardabweichung statistisch kompatibel mit Null ($B_{\text{offset}} \approx 0\ \text{mT}$). Dies spricht für eine hervorragende Kompensation externer Störfelder oder eine hochwirksame magnetische Abschirmung (z.B. mittels einer $\mu$-Metall-Hülle) während der Messung.

3. Periodizität $\Delta B$ und effektive Junction-Fläche $A_{\text{eff}}$

Die Periodizität der Oszillation, welche dem Durchtritt eines magnetischen Flussquants $\Phi_0$ durch die effektive Barrierefläche entspricht, beträgt $\Delta B = (1{,}0907 \pm 0{,}0191)\ \text{mT}$.

Unter Verwendung des magnetischen Flussquants $\Phi_0 = \frac{h}{2e} \approx 2{,}0678 \times 10^{-15}\ \text{Vs}$ lässt sich die effektive Fläche $A_{\text{eff}}$ des Josephson-Kontakts über die Beziehung

$$A_{\text{eff}} = \frac{\Phi_0}{\Delta B}$$

bestimmen. Daraus ergibt sich eine effektive Fläche von:

$$A_{\text{eff}} = (1{,}8959 \pm 0{,}0332)\ \mu\text{m}^2$$

Diese effektive Fläche $A_{\text{eff}}$ setzt sich geometrisch aus der physikalischen Barrierebreite $w = 10\ \mu\text{m}$, der Barrierendicke $D = 30\ \text{nm}$ und den London'schen Eindringtiefen $\lambda_{L,1}$ und $\lambda_{L,2}$ der beiden supraleitenden Elektroden zusammen:

$$A_{\text{eff}} = w \cdot (D + \lambda_{L,1} + \lambda_{L,2})$$

Da die reine geometrische Barriere $D$ im Bereich weniger Nanometer liegt, wird die effektive Fläche dominant durch das Eindringen des Magnetfeldes in die supraleitenden Randbereiche bestimmt.

4. London'sche Eindringtiefe $\lambda_L$

Unter der physikalisch gerechtfertigten Annahme symmetrischer supraleitender Elektroden ($\lambda_{L,1} = \lambda_{L,2} = \lambda_L$) lässt sich die London'sche Eindringtiefe direkt aus der effektiven Fläche $A_{\text{eff}}$ isolieren:

$$\lambda_L = \frac{1}{2} \left( \frac{A_{\text{eff}}}{w} - D \right)$$

Mit den Probenabmessungen von $w = 10\ \mu\text{m}$ und $D = 30\ \text{nm}$ ergibt sich ein Wert von:

$$\lambda_L = (79{,}79 \pm 1{,}66)\ \text{nm}$$

Die London'sche Eindringtiefe beschreibt die charakteristische Skala, auf der ein äußeres Magnetfeld exponentiell im supraleitenden Material abgeschirmt wird ($B(z) = B_0 e^{-z/\lambda_L}$). Der ermittelte Wert von ca. $80\ \text{nm}$ steht in hervorragender quantitativer Übereinstimmung mit Literaturwerten für dünne Niob-Schichten (Nb), welche typischerweise zwischen $80\ \text{nm}$ und $110\ \text{nm}$ liegen. Geringe Abweichungen zu massiven Nb-Einkristallen ($\lambda_L \approx 40\ \text{nm}$) lassen sich physikalisch durch eine reduzierte mittlere freie Weglänge der Elektronen in gesputterten Dünnfilmen (Dirty-Limit) erklären.

5. Kritische Stromdichte $j_c$ und Josephson-Eindringtiefe $\lambda_J$

Die kritische Stromdichte $j_c$ normiert den maximalen Josephson-Strom auf die rein geometrische Kontaktfläche ($A_{\text{geo}} = w \cdot L = 100\ \mu\text{m}^2$):

$$j_c = \frac{I_0}{w \cdot L} = (5909{,}30 \pm 352{,}51)\ \text{A/cm}^2$$

Diese hohe Stromdichte ist charakteristisch für qualitativ hochwertige, dünne Tunnelbarrieren. Aus $j_c$ und der magnetischen Dicke $d_{\text{mag}} = D + 2\lambda_L$ lässt sich die Josephson-Eindringtiefe $\lambda_J$ berechnen:

$$\lambda_J = \sqrt{\frac{\Phi_0}{2\pi \mu_0 d_{\text{mag}} j_c}} = (48{,}35 \pm 1{,}50)\ \mu\text{m}$$

Die Josephson-Eindringtiefe stellt die fundamentale Längenskala dar, auf der sich der supraleitende Phasenunterschied $\varphi(x)$ entlang der Barriere räumlich anpassen kann. Gleichzeitig beschreibt sie die Abschirmungstiefe des durch den Josephson-Strom selbst induzierten Magnetfeldes (Selbstabschirmung).

6. Klassifizierung des Kontakts: Das Verhältnis $L/\lambda_J$

Ein zentrales Kriterium für die Interpretation der Interferenzmessungen ist das Verhältnis der Kontaktlänge $L$ zur Josephson-Eindringtiefe $\lambda_J$:

$$\frac{L}{\lambda_J} = 0{,}2068 \pm 0{,}0064$$

Da gilt:

$$\frac{L}{\lambda_J} \approx 0{,}21 \ll 1$$

befindet sich der untersuchte Josephson-Kontakt im sogenannten "Short-Junction-Limit" (kurzer Kontakt). Diese Einordnung hat weitreichende physikalische Konsequenzen für die Validität unserer Auswertung:

Homogene Stromverteilung: Da die geometrische Ausdehnung $L = 10\ \mu\text{m}$ deutlich kleiner als die charakteristische Abschirmungslänge $\lambda_J \approx 48\ \mu\text{m}$ ist, kann der fließende Josephson-Strom das im Kontakt herrschende Magnetfeld nicht nennenswert abschirmen (vernachlässigbare Selbstabschirmungseffekte).

Linearer Phasenverlauf: Der räumliche Verlauf der supraleitenden Phase $\varphi(x)$ entlang des Kontakts ist in exzellenter Näherung linear und wird ausschließlich durch das homogene externe Magnetfeld diktiert ($\frac{\partial \varphi}{\partial x} = \frac{2\pi d_{\text{mag}}}{\Phi_0} B_{\text{ext}}$).

Validität des Fraunhofer-Modells: Diese Linearität ist die fundamentale mathematische Voraussetzung für die mathematische Fourier-Transformation, die zur idealen Fraunhofer-Form ($I_c \propto \vert{}\text{sinc}(x)\vert{}$) führt. Wäre $L \ge \lambda_J$, würden nichtlineare Effekte (wie die Entstehung von Josephson-Wirbeln oder Flussquanten-Eindringen) das Interferenzmuster stark verzerren, was sich in asymmetrischen, hysteretischen Kurven äußern würde. Die exzellente Anpassbarkeit unserer Daten bestätigt das Vorliegen des Short-Junction-Limits experimentell.

7. Konstanter Strom-Offset $C$

Der Fit liefert einen kleinen, aber statistisch signifikanten konstanten Strom-Offset von $C = (9{,}29 \pm 1{,}05)\ \mu\text{A}$ (bzw. $0{,}0093 \pm 0{,}0010\ \text{mA}$). Mathematisch verhindert dieser Parameter, dass die Minima (Nullstellen) des Fraunhofer-Musters exakt auf die Nulllinie absinken. Physikalisch lässt sich dieses Verhalten wie folgt erklären:

Fluss-Inhomogenitäten: Eine nicht perfekt homogene Stromverteilung im Übergang führt dazu, dass sich die Phasenunterschiede über die Breite des Kontakts nicht in allen Punkten exakt destruktiv interferieren.

Thermische Rauscheffekte: Thermische Fluktuationen (Johnson-Rauschen) im Messkreis können ein verfrühtes Zurückspringen aus dem supraleitenden Zustand in den resistive Zustand induzieren.

Leckströme: Ein kleiner normalleitender Parallelpfad (Shunt) über den Rand der Barriere trägt zu einem ohmschen Hintergrundstrom bei.

Unsicherheitsbetrachtung und kritische Diskussion

Die statistischen Unsicherheiten der Fitparameter wurden über die Kovarianzmatrix der nichtlinearen Ausgleichsrechnung bestimmt. Sie reflektieren die Streuung der Messpunkte um die theoretische Kurve.

Als wesentliche systematische Fehlerquellen sind zu nennen:

Kalibrierung des Spulenpaares: Der Umrechnungsfaktor von Spulenspannung zu Magnetfeldstärke ($1{,}927\ \text{V} \rightarrow 5{,}8\ \text{mT}$) unterliegt Fertigungstoleranzen der Helmholtz-Spulen (Geometrieabweichungen, Windungsdichte). Eine systematische Abweichung hierbei verschiebt die Periodizität $\Delta B$ direkt proportional und führt zu einem systematischen Fehler in der Berechnung von $A_{\text{eff}}$ sowie nachfolgend $\lambda_L$ und $\lambda_J$.

Kompensationsströme und thermische Spannungen: An den Kontakten der Zuleitungen im Kryostaten entstehen durch Temperaturgradienten thermoelektrische Spannungen (Seebeck-Effekt). Diese verschieben die gemessenen Spannungen geringfügig, was sich in der Auswertung durch den Plateau-Offset im Bereich von ca. $-0{,}09\ \text{mV}$ äußert. Da dieses Offset jedoch bei der Differenzbildung der kritischen Ströme rechnerisch eliminiert wurde, ist der Einfluss auf $I_0$ minimiert.

Phasenrauschen und externe HF-Einkopplung: Hochfrequente elektromagnetische Störungen aus der Laborumgebung können wie ein effektives magnetisches Rauschen wirken, welches die Minima des Oszillationsmusters systematisch verschmiert und somit künstlich den Offset-Parameter $C$ erhöht.

Zusammenfassend zeigt die hervorragende Übereinstimmung der Messdaten mit dem theoretischen Fraunhofer-Modell (ersichtlich an den geringen relativen Fehlern der Hauptparameter $I_0$ und $\Delta B$ von unter $6\ \%$), dass der Josephson-Kontakt im untersuchten Feldbereich makroskopisch quantenmechanisch kohärent arbeitet und die fundamentale Flussquantisierung im Supraleiter erfolgreich demonstriert wurde.