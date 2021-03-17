# nsoFI
***********************************************
NSO.py  Night Sky Objects - Yötaivaan kohteet

Python 3.7: ssä tehty ohjelma tähtien päivittäisen näkyvyyden tarkistamiseksi Kanariansaarten ja Chilen Slooh-observatorioissa sekä eräissä Suomen kaupungeissa. Ohjelma ei tarkista päivän valoisuutta.

Vaadittavat tiedostot:
================
1. NSO.py-ohjelma
2. - 4. MessierObjects.csv, CaldwellObjects.csv ja VariableStarObjects.csv: pilkuin erotetut arvot - tiedosto taivaskoordinaateista (rektakensio ja deklinaatio). Voit lisätä omia mukautettuja tietoja tiedostoon, jos haluat
5. Observatory.csv: observatorion koordinaattien tiedosto (nimi, pituus, leveys, aloitus- ja päättymisajat)
6. Mallikuva.jpg: aloitusnäyttö, voidaan jättää pois, et vain näe kuvaajaa ohjelman käynnistyessä

Ohjelman suorittaminen:
====================
Käynnistä NSO.py python-alustallasi
- valitse kohde ja observatorio
- Jos haluat tallentaa kuvaajan, muuta valintaruudun tilaa. Tiedostonimessä on esineen nimi, observatorion nimi ja päivämäärä.
- Ohjelmarakenteen vuoksi ei tallenna -vaihtoehto tallentaa kuvaajan muodossa 'sivutallennus.jpg'

Python-moduulit:
===============
Kaikki Python-moduulit ovat yleensä Python 3.7 -ympäristössä, mutta ne on ehkä asennettava.
- panda
- matplotlib
- pyqt5
- math
- datetime
- re
- csv

Versiohistoria:
================
Ver 1.2 10.06.2020
- lisäsi Messier, Caldwell ja Variable Stars -valintalaatikon

Ver 1.1 - 03.06.2020
- lisäsi kuvan tallennusvaihtoehdon
- korjattu päivämäärän selitys x-akselilla

Ver 1.0 30.5.2020
- ensimmäinen versio


Kiitokset:
==========
Tähtitieteelliset laskelmat perustuvat Paul Schlyterin verkkosivustoon https://stjarnhimlen.se/comp/ppcomp.html
Inspiraatio Jarmo Ruuthin verkkosivustolta https://ruuth.xyz/AstroMosaic.html
