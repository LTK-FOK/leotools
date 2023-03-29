Ez a Lechner Earth Observation Toolset (leotools), egy python csomag, amelyet a Lechner Tudásközpont Űrtávérzékelési Osztálya tart fent és használ.

romeo.ife.komolafe@lechnerkozpont.hu

# Telepítés

Telepítsd a leotools csomagot egy Python környezetbe conda vagy pip csomagkezelők segítségével!

## conda

Telepítéshez a letöltött vagy az előzőleg buildelt `.tar.bz2` állományt adjuk meg.

```bash
conda install <dist>
```

Ha buildelni szeretnénk, először klónozzuk a csomag repository-ját. A conda környezetnek tartalmaznia kell a `conda-build`, `setuptools` és `wheel` csomagokat.

```bash
git clone https://github.com/LTK-FOK/leotools
conda build leotools
```

## pip

Telepíthetünk közvetlenül GitHubról.

```bash
pip install git+https://github.com/LTK-FOK/leotools
```

Telepíthetjük a letöltött vagy az előzőleg buildelt `.whl` vagy `tar.gz` állományt is.

```bash
pip install <dist>
```

Telepíthetünk közvetlenül a klónozott repository-ból.

```bash
git clone https://github.com/LTK-FOK/leotools
cd leotools
python setup.py install
```

Ha buildelni szeretnénk, a Python környezetnek tartalmaznia kell a `setuptools` és `wheel` csomagokat.

```bash
git clone https://github.com/LTK-FOK/leotools
cd leotools
python -m build
```

# Használat

A leotools több módon is használható:

1. Egyes almodulokat és eszközöket a megszokott módon lehet más Python szkriptekbe importálni.
2. Bizonyos függvények a parancssorból is elérhetőek.

## Parancssoros interfész

Ha a leotools telepítési környezete aktív, bizonyos függvények elérhetőek a parancssorból.

Megjeleníthetjük a súgót a `help` paranccsal.

```bash
leotools help
```

Az elérhető függvények listázhatók a `list` paranccsal.

```bash
leotools list
```

Az egyes függvények dokumentációja is kiiratható.

```bash
leotools help preproc
```

Hasznos lehet kiiratni általános információkat gyakran használt argokról.

```bash
leotools info
leotools help info
```

Egy függvény így használható. Opcionális argokat kötőjelek segítségével adhatunk meg.

```bash
leotools tile example.zip outdir -meta_dir outdir
```

## Importálási konvenciók

Kiirathatjuk a csomag különböző információit.

```python
import leotools
print(leotools.__description__)
# Lechner Earth Observation Toolset, containing crucial classes and functions.
```

Az elérhető eszközöket a tartalmazó modulból kell importálnunk, azok nem importálhatók közvetlenül a csomagból.

```python
from leotools import gistools as ltgt
from leotools.basetools import ProcessTimer, load_files
```

Egy függvény vagy osztály belső dokumentációja a `__doc__` tulajdonságban van tárolva. Ez listázza az összes paramétert, amit megadhatunk neki.

```python
print(image_to_array.__doc__)
# Loads selected bands of an image into a numpy array.
#
# The dimensions of the output array are `count, rows(y), cols(x)`, or
#     `rows(y), cols(x)` for a single-band image.
#
# Args:
#     image (path): The image to load.
#     bands (int or list, optional): Selected bands. Bands can be omitted,
#         duplicated and rearranged.
#
# Returns
#     array (ndarray): The loaded array.
```

# Segédfunkciók

A kiegészítő eszközök, mint a fájlkezelés és az időmérés a `basetools` modulban vannak összegyűjtve.

## Útvonalak ellenőrzése

Ellenőrizhetjük, hogy egy útvonal létezik-e. Ha létezik, `mode`-tól függetlenül semmi nem történik.

```python
path = "path/to/example/file.tif"
check_path(path)
### Semmi nem történik
```

Nem létező útvonal esetén a `mode` arg határozza meg a működést.

```python
path = "path/to/nonexistent/file.tif"
check_path(path, mode=1)
# FileNotFoundError: This path does not exist: path/to/example/folder
```

|  mode  | Viselkedés                                                   |
| ------ | ------------------------------------------------------------ |
| 0      | Hibát dob, ha az útvonal nem létezik.                        |
| 1      | Létrehozza az útvonal végpontját, de hibát dob, ha a szülők hiányoznak. |
| 2      | Létrehozza az útvonal összes hiányzó elemét.                 |

## Útvonalak betöltése

Útvonalakat több forrásból tölthetünk be a `load_files` függvény segítségével. Érdemes abszolút útvonalakat használni. Ha az útvonalban visszaperek vannak, használjunk `r"D:\files\subdir\ex4.tif"` formátumú "raw" stringet.

Nézzük a következő mappastruktúrát:

```
D:/files
  ├─ ex1.tif
  ├─ ex2.jp2
  └─ subdir
      ├─ ex3.tif
      └─ ex4.tif
paths.txt
```

Megadhatunk egy útvonalat vagy útvonalak listáját.

```python
single_path = load_files("D:/files/subdir/ex4.tif")
print(single_path)
# [WindowsPath('D:/files/subdir/ex4.tif')]

multiple_paths = load_files(["D:/files/ex1.tif", "D:/files/subdir/ex3.tif"])
print(multiple_paths)
# [WindowsPath('D:/files/ex1.tif'), WindowsPath('D:/files/subdir/ex3.tif')]
```

Ha egy mappát adunk meg, annak tartalmát listázza ki.

```python
dir_contents = load_files("D:/files/subdir")
print(dir_contents)
# [WindowsPath('D:/files/subdir/ex3.tif'), WindowsPath('D:/files/subdir/ex4.tif')]
```

Ha egy .txt fájlt adunk meg, minden sort külön útvonalként kezel. Az üres sorokat és a kettőskereszttel kezdett kommenteket ignorálja.

A `paths.txt` fájl tartalma:

```
D:/files/ex1.tif
D:/files/ex2.jp2 # Ez egy komment

D:/files/subdir/ex3.tif
```

```python
txt_contents = load_files("D:/paths.txt")
print(dir_contents)
# [WindowsPath('D:/files/ex1.tif'), WindowsPath('D:/files/ex2.jp2'), WindowsPath('D:/files/subdir/ex3.tif')]
```

Szűrhetünk bizonyos fájlokra a `pattern` változtatásával. A `recursive` argot bekapcsolva a függvény az almappákat is átfésüli.

```python
paths = load_files("D:/files", pattern='*.tif', recursive=True)
print(dir_contents)
# [WindowsPath('D:/files/ex1.tif'), WindowsPath('D:/files/subdir/ex3.tif'), WindowsPath('D:/files/subdir/ex4.tif')]
```

Az `str_out` argot `True`-ra állítva az alapértelmezett `pathlib.Path` objektum helyett az útvonalakat stringként kaphatjuk vissza.

## Időzítés

Időzítésre a `ProcessTimer` osztály szolgál. Amikor elindítjuk vagy megállítjuk az időzítőt, megadhatunk neki egy nevet, amit kiír, amikor megáll.

```python
from time import sleep

timer = ProcessTimer()
timer.start("My time")
sleep(1)

timer.pause()
sleep(1) ### Itt megállítottuk, ez nem fog beleszámítani a végső időbe

timer.start()
sleep(1)
timer.stop()
# My time: 2.00s
```

Az időzítőket egymásba ágyazhatjuk. Ha elindítunk egy már futó időzítőt, az automatikusan megállítja magát, mielőtt újraindul.

```python
timer = ProcessTimer()
subtimer = ProcessTimer()

timer.start("Full process")
for i in range(3):
    subtimer.start()
    myfunc(i)

subtimer.stop()
timer.stop()

# 10.12s
# 9.89s
# 9.56s
# Full process: 29.57s
```

Ha nem szeretnénk kiírni az időt, csak változóban visszakapni, állítsuk `False`-ra a `verbose` argot.

```python
timer = ProcessTimer(verbose=False)
timer.start("My time")
myfunc()
process_time = timer.stop() ### Nem írja ki

print(process_time)
# My time: 1m 42.56s
```

# Térinformatikai operációk

Az általános térinformatikai eszközök a `gistools` modulban vannak összegyűjtve.

## Profil és kép információk

Egy kép profilja kinyerhető dictionary formában. Ez vetületi, tömörítési és más nem spektrális információkat tartalmaz. Egyéb adatokat -- metaadat címkék, fájl méret, olvasási idő -- is megkaphatunk a megfelelő függvényekkel.

```python
profile = get_profile("example.tif")
print(profile)
# {'driver': 'GTiff', 'dtype': 'uint16', 'nodata': 0.0, 'width': 7990, 'height': 8100, 'count': 7, 'crs': CRS.from_epsg(23700), 'transform': Affine(30.0, 0.0, 776100.0,
#       0.0, -30.0, 517800.0), 'blockxsize': 256, 'blockysize': 256, 'tiled': True, 'compress': 'zstd', 'interleave': 'band'}

tags = get_tags("example.tif")
print(tags)
# {'AREA_OR_POINT': 'Area', 'comp_level': '1', 'comp_pred': '2', 'date': '20220311', 'path': '186', 'refl': 'boa', 'sat': 'L8', 'sources': 'LC08_L2SP_186026_20220311_20220321_02_T1'}

file_size = get_file_size("example.tif")
print(file_size)
# 0.43994333129376173

read_time = get_read_time("example.tif")
print(read_time)
# 3.14s
```

Ezek az adatok rögtön kiirathatók formázva a `print_profile` függvénnyel. Az olvasási idő lemérése eltarthat egy ideig, ezt kikapcsolhatjuk a `read_time` arg `False`-ra állításával.

```python
print_profile("example.tif", read_time=False)
# example.tif
#
# driver:  GTiff
# dtype:  uint16
### ...
# sat:  L8
# sources:  LC08_L2SP_186026_20220311_20220321_02_T1
#
# File size: 0.440 GB
```

Átállíthatjuk egy kép nodata értékét a `set_nodata` függvénnyel.

```python
set_nodata("example.tif", 999) ### 999-re állítjuk a nodata értéket
```

Egy kép vetületi információit a `strip_proj` függvénnyel törölhetjük.

```python
print(get_profile("example.tif")['crs'])
# EPSG:23700
strip_proj("example.tif")
print(get_profile("example.tif")['crs'])
# None
```

Egy kép metaadat címkéit a `strip_tags` függvénnyel törölhetjük.

```python
print(get_tags("example.tif"))
# {'AREA_OR_POINT': 'Area', 'comp_level': '1', 'comp_pred': '2', 'date': '20220311', 'path': '186', 'refl': 'boa', 'sat': 'L8', 'sources': 'LC08_L2SP_186026_20220311_20220321_02_T1'}
strip_tags("example.tif")
print(get_tags("example.tif"))
# {'AREA_OR_POINT': 'Area'}
```

## Kép és tömb

Betölthetünk egy képet egy numpy tömbbe az `image_to_array` függvénnyel. A `bands` arggal megadhatjuk a betöltött sávokat. Ha nem adjuk meg, a kép összes sávja betöltésre kerül eredeti sorrendben. A sávokat bármilyen sorrendben, akár duplikálva is megadhatjuk.

```python
image = "example.tif"

array1 = image_to_array(image) ### Az összes sáv
print(array1.shape)
# (7, 8100, 7990)

array2 = image_to_array(image, bands=[1, 2, 4, 3, 4]) ### Tetszőleges sávkombináció
print(array2.shape)
# (5, 8100, 7990)

array3 = image_to_array(image, 1)
print(array3.shape)
# (8100, 7990)

array4 = image_to_array(image, [1])
print(array4.shape)
# (1, 8100, 7990)
```

Egy numpy tömböt menthetünk képként az `array_to_image` függvénnyel. A mentendő kép profilja a `sample_image` és a `profile` argokból áll elő. Legalább az egyikre szükség van, utóbbi felülírja az előbbit.

```python
profile = {'compress': 'lzw'} ### Megváltoztatjuk a tömörítési algoritumst
array_to_image(array, "output.tif", sample_image="example.tif", profile=profile)
### Eredmény: output.tif
```

## Átkódolás és extent kerekítés

Egy tömböt átkódolhatunk, kicserélhetjük bizonyos értékeit. Megadhatjuk a `default_value` argot, ha a nem megadott értékeket kicserélnénk egy bizonyos értékre.

```python
array = np.array([
    [1, 2, 3],
    [1, 2, 3],
    [4, 5, 6]
])

new_array = remap_array(array, {1: 10, 2: 20, 3: 30})
print(new_array)
# [[10 20 30]
#  [10 20 30]
#  [ 4  5  6]]

new_array = remap_array(array, {1: 10, 2: 20, 3: 30}, default_value=0)
print(new_array)
# [[10 20 30]
#  [10 20 30]
#  [ 0  0  0]]
```

Egy kép extentjeinek listáját tetszőleges felbontású gridre kerekíthetjük.

```python
old_extent = [12, 200, 315, 1189] ### xmin, ymin, xmax, ymax
new_extent = round_extent(old_extent, 100)
print(new_extent)
# [0, 200, 400, 1200]
```

## Maszkolás

Egy tömböt kimaszkolhatunk egy másik tömbbel a `mask_array` függvény segítségével. Alapvetően az `1` vagy `True` értékkel rendelkező cellákat tartja meg a függvány.

```python
array = np.array([
    [1, 2],
    [3, 4]
])

mask = np.array([
    [1, 0],
    [0, 1]
])

masked_array = mask_array(array, mask, nodata=0)
print(masked_array)
# [[1 0]
#  [0 4]]
```

Egy raszeteres képet kimaszkolhatunk vektoros állományokkal a `multi_mask` függvény segítségével. Utóbbi állományok lehetnek fájlok és geometriák is. Alapesetben a geometriák alatti területet tartja mega függvény.

```python
import geopandas as gpd

area1 = gpd.read_file("area1.shp") ### geopandas GeoSeries objektum

multi_mask(
    "example.tif",
    "masked.tif",
    nodata=0,
    area_masks=[area1, "area2.shp"], ### A terület maszkok
    bound_masks=["bound.shp"] ### A körvonal maszkok
)
### Eredmény: masked.tif
```

Mindkét függvény eredményét megfordíthatjuk az `invert` arggal.

```python
array = np.array([
    [1, 2],
    [3, 4]
])

mask = np.array([
    [1, 0],
    [0, 1]
])

masked_array = mask_array(array, mask, nodata=0, invert=True)
print(masked_array)
# [[0 2]
#  [3 0]]
```

## Sávok kinyerése és indexek előállítása

Kinyerhetünk bizonyos sávokat egy képből az `extract_bands` fügvénnyel. Ezek a sávok külön egysávos képként jelennek meg a kimeneti mappában.

```python
extract_bands("example.tif", "output_dir", bands=[1, 3, 5])

### Eredmények:
### output_dir/example_1.tif
### output_dir/example_3.tif
### output_dir/example_5.tif
```

Előállíthatunk indexeket, ha megadjuk a nevüket és a képletüket egy dictionary-ben. A képletekben használhatjuk az `ndiff` függvényt normalizált különbség értékek egyszerű előállításához.

```python
bands = {
    'b3': "b[3]", ### A 3-as sáv
    'ndvi': r"(b[8]-b[4]) / (b[8]+b[4])", ### Index
    'ndvi_alt': r"ndiff(b[8], b[4])" ### Az `ndiff` függvény használata
}

extract_bands("example.tif", "output_dir", bands=bands, dtype='uint16')

### Eredmények:
### output_dir/example_b3.tif
### output_dir/example_ndvi.tif
### output_dir/example_ndvi_alt.tif
```

Ha a számított indexek adattípusa más lenne, mint a megadott `dtype`, rakjunk egy átszámítást a képletbe.

```python
bands = {'ndvi': r"(ndiff(b[8], b[4])+1)*10000"}
extract_bands("example.tif", "output_dir", bands=bands, dtype='uint16')
```

A `constants` modulban megtalálhatók a Landsat és Sentinel-2 felvételek gyakran használt indexeinek nevei és képletei dictionary formában.

```python
from leotools.constants import NDXI_L8, NDXI_S2

print(NDXI_L8['tct1'])
# b[1]*0.3029 + b[2]*0.2786 + b[3]*0.4733 + b[4]*0.5599 + b[5]*0.5080 + b[6]*0.1872

print(NDXI_S2['grvi'])
# b[5] / b[3]
```

## Stackelés

Több képet stackelhetünk a `stack_images` függvénnyel.

```python
stack_images("input_dir", "stacked.tif")
print(get_profile("stacked.tif")['count'])
# 21
```

A `band` arggal megadhatjuk az input képeknek egy bizonyos sávját, amit össze akarunk rakni. Ha megadjuk a `bands_csv` argot, kiirathatunk egy csv fájlt, ami listázza az egyes sávok forrását.

```python
stack_images("input_dir", "stacked.tif", band=4, bands_csv="bands.csv")
print(get_profile("stacked.tif")['count'])
# 3
```

A bands.csv tartalma:

```
1,D:\input_dir\ex1.tif,4
2,D:\input_dir\ex2.tif,4
3,D:\input_dir\ex3.tif,4
```

Ha tömböket szeretnénk stackelni, használjuk a `stack_arrays` függvényt. Ennek eredményét elmenthetjük képként az `array_to_image` függvény segítségével. Használhatunk egyszerre 2 dimenziós (y, x) és 3 dimenziós (sávszám, y, x) tömböket is.

```python
array1 = image_to_array(r"D:\romeo\out\20220311_L8_T1_186_boa_eov.tif", 1)
print(array1.shape)
# (8100, 7990)

array2 = image_to_array(r"D:\romeo\out\20220311_L8_T1_186_boa_eov.tif", [2])
print(array2.shape)
# (1, 8100, 7990)

array3 = image_to_array(r"D:\romeo\out\20220311_L8_T1_186_boa_eov.tif", [3, 4, 5])
print(array3.shape)
# (3, 8100, 7990)

stacked_array = stack_arrays([array1, array2, array3])
print(stacked_array.shape)
# (5, 8100, 7990)
```

## OVR és AUX fájlok előállítása

Gyárthatunk külső overview-t, ami a piramisrétegeket tartalmazza.

```python
make_ovr("example.tif", levels=[2, 4, 8])
### Eredmény: example.tif.ovr
```

A kép spektraális statisztikáit kiírhatjuk egy .aux fájlba.

```python
make_aux("example.tif")
### Eredmény: example.tif.aux.xml
```

# Előfeldolgozás

A műholdfelvételek előfeldolgozásához használt eszközök a `preproc` modulban vannak összegyűjtve.

## Általános esetben

Az előfeldolgozás minden lépését kezeli a `preproc` függvény. Ez a .zip és .tar archívumokból mozaikolt .tif datatakeket állít elő. Megadhatjuk a `meta_output_dir` argot, és a függvény ide fogja átrakni a metaadatokat tartalmazó szöveges állományokat. Megadhatjuk a `temp_dir` argot, ellenkező esetben a függvény maga kezeli az ideiglenes fájlok elhelyezését és törlését.

```python
preproc("input_dir", "output_dir", meta_output_dir="meta_dir")
```

### Tile-ok feldolgozása

Ha az archívumokból képeket szeretnénk csinálni mozaikolás nélkül, használjuk a `reproj_tile` függvényt.

Az `ls_kwargs` és `s2_kwargs` argok segítésgével argokat továbbíthatunk a Landsat és Sentinel-2 feldolgozásoknak. A `used_bands` arggal meghatározhatjuk, hogy melyik sávokat szeretnénk látni a feldolgozott képben. Az `incl_gt` arg Sentinel-2 felvételeknél a generálás időpontját a név mögé rakja.

```python
l8_zip = "LC08_L2SP_186026_20220311_20220321_02_T1.tar"
s2_zip = "S2A_MSIL2A_20220102T095411_N0301_R079_T34TCR_20220102T114458.zip"

reproj_tile(l8_zip, "output_dir")
### Eredmény: output_dir/20220311_L8_T1_186_26_boa_eov.tif

### Ha ssak az első 3 sávot használjuk, a 
reproj_tile(
    [l8_zip, s2_zip],
    "output_dir",
    ls_kwargs={'used_bands': ['B1', 'B2', 'B3']},
    s2_kwargs={'used_bands': ['B01', 'B02', 'B03']}
)
### Eredmények:
### output_dir/20220311_L8_T1_186_26_boa_eov.tif
### output_dir/20220102_s2a_r079_tcr_boa_eov.tif


reproj_tile(s2_zip, "output_dir", s2_kwargs={'incl_gt': True})
### Eredmény: output_dir/20220102_s2a_r079_tcr_boa_eov_20220102T114458.tif
```

### Feldolgozott tile-ok mozaikolása

Ha már meglévő képeket szeretnénk mozaikolni, használjuk a `merge_datatake` függvényt. Fonotos, hogy ez csak a `reproj_tile`-al készült képeken tud dolgozni.

```python
merge_datatake("input_dir", "output_dir", meta_dir="input_dir", meta_output_dir="output_dir")
### Itt ugyanazokat a mappákat használjuk a képekhez és a metadata fájlokhoz
```

## Meglévő képek frissítése

Egy képet a leotools által használt egységes forma szerint átalakíthatunk a `reformat` függvénnyel. A kép spektrális és vetületi értékei nem változnak. A `compatible` arg bekapcsolásával egy kevésbé gazdaságos, de széles körben kompatibilis formátumot használ a függvény.

```python
reformat("example.tif", "output_dir", compatible=False, suffix='')
### Eredmény: output_dir/example.tif
```

A `make_extras` függvény egyesíti a `make_ovr` és `make_aux` fájlok funkcionalitását.

```python
make_extras(input_paths, ovr=True, aux=True, recursive=False)

### Eredmények:
### example.tif.ovr
### example.tif.aux.xml
```

## Általános paraméterek

A `preproc` modul sok függvényében elérhetőek bizonyos univerzális argok. A `recursive` arg meghatározza, hogy a bemeneti mappák rekurzívan vannak-e átfésülve a [`load_files`](#útvonalak-betöltése) függvény szerint. A `mode` arg a kimeneti mappákat ellenőrzi a [`check_path`](#útvonalak-ellenőrzése) függvény szerint.

## Képek szűrése címkék alapján

Az elkészült képek listáját megszűrhetjük a metaadat címkéik alapján.

A `list_uniques` függvény segítségével kilistázhatjuk egy-egy címke egyedi értékeit.

```python
file_list = load_files(r"D:\image_dir")

uniqe_path = list_uniques(file_list, 'path')
print(unique_path)
# ['r036' 'r079' 'r136']
```

A `filter_images` segítségével szűrhetünk képekre. Az `expressions` argban adhatjuk meg a szűrő kifejezések listáját `"<tag> <operator> <value>"` formában pl. `"date >= 20220101"`. A szokásos `==`, `!=`, `>`, `<`, `>=` és `<=` operátorok mellett a `cc` és `!c` operátorok is elérhetőek, ezek a "tartalmaz" és "nem tartalmaz" operációkat jelentik.

```python
filtered_list = filter_images(file_list, ["date >= 20210201", "date < 20210301", "path cc 36"])
print(filtered_list)
# [WindowsPath('D:/image_dir/20220203_s2b_r036_boa_eov.tif')
#  WindowsPath('D:/image_dir/20220208_s2a_r036_boa_eov.tif')
#  WindowsPath('D:/image_dir/20220213_s2b_r036_boa_eov.tif')
#  WindowsPath('D:/image_dir/20220218_s2a_r036_boa_eov.tif')
#  WindowsPath('D:/image_dir/20220210_s2b_r136_boa_eov.tif')
#  WindowsPath('D:/image_dir/20220225_s2a_r136_boa_eov.tif')]
```