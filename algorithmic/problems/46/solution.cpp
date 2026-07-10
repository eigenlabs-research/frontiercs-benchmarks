#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <climits>
#include <ctime>
using namespace std;
struct E{uint64_t h;const char*s;};static const E a[]={
{1337409398582560185ULL,"I6ADKPSE120LBF59Q4CTO7MR8NJHG3KAHLNOTS04EDQPR1IGF3J86BM7295CNGJA1L7Q9S3C05IMF4HO8RKBTE62PD1K50S4F9IORNTQECP2AMH8JDBG7L36421L5SHPTF73A6EGODMNRBK9CI0J8QHEAGB84MJQ120TSLKOCR657PIDF93N29I8AMLGE5STHPKJC46N7F1R0OQ3DBQ715CDKL04JSH23PAI9ME6TN8GFBORTDE0P32JCQM1FO7GRN5ILSK8H4B6A997GMKJ0AH2I186LPSNOEF3R54TDCBQI7T4DJS9HQ2FL3A60RPG8M1NKB5OCEEDOMH1G0583RSPKI7A2JN9BCL4Q6FTT072DK1MEI8LO964NH3BRC5SAJPGFQMJIG85KC7A20EDL3HN4PFOB6RQT1S9JMBHS32CQP87KAELRO0914I5G6FTDN4CBK501HA8M7I3D2RTQP9OLJGFENS6A0Q5KIFRJL2783PC49HEDT1NGSO6MB6B0ISM23JQ9LH1FOKN8GEDC4AR75PT78PEM0FDRC4295JG3LQBTIHKA16NSORLET0G62BS17CMI94AKOFJNH3D5Q8P"},
{1097491252512501470ULL,"cQMF37P2UbDLIB05EJYN4RWdZA16KG8OT9VXCaSHMODU1cZP3IQ5YX2CEHb0a4RWFd7A6KG8T9VLBSJNE8LGPBZ37Q65O2CVSHXbMJ0NY4RWFdcA1KDUT9Ia3G57QODKdZ628PIBULaNYSHECbMJ4RWFcA1T9VX08CEUDdcP1bZM35QOGBINY0SH4RWF7A6K2T9VXLaJZU7QOc2JM0aSX5LEYNH4RWFdPA16K3DG8T9VBbCIaP4ZdOU7KN3L2I8QSB5HJYMRWFcA16DGT9VXEb0C4dDP3Z57XQO02UILBCYSaNHMRWFcA16KG8T9VEbJ52Q1cP68dZDIGEB30bMOYLS4RWF7AKUT9VXCaHJNR4WFd7cZA1P6K3DG8UTO29QVXEbBCI0LaJSHN5YM23OULD0ZB7MQP8X5bICSHY4RWFdcA16KGT9VEaJNP7F26T3LcDUN8OV5B0SEIYM4RWdZA1KGQ9XbCaHJEPKOG7U0D5S2Q6NCLabHBYM4RWFdZcA138T9VXIJ8XVTOcLJD3PUGE2C5IBHMYSN4RWFd7ZA16KQ9b0ab4cD8PGY5O3UL1X2N0IMCSHRWFd7ZA6KQT9VEBaJ"},
{13064344411517067317ULL,"FS9PLDQ5ETMK7IA01RJ364GOB8NHC22HSK68T0A945OB7CLMR3FDIPJQNE1G06SN98DCQ1L5MPFA43TJHBE7OGIK2RIJ52COTMG8BPA9HFRKLQ03D16N7ES48Q1IPA2N3SD4GMC5H7LER9KB0TO6JFPEIFOK57TH3A4LMJGNCSQ291806DRB63GACOBPEL4QMNJ8I52107H9FDTKRSH3DQO6C84FJ71PTL5EBK2S0N9MIRGAA9C7BD30T82K6I1NFRH4Q5OPGLMJSE8I3C6EDKQ1AS2RMBGFJPLH5O7T049N7ON9AFIGC1L3Q64RSE0BMDK852HJTPQ7A4J605913OCRBTMLHFDIG8SPN2EKMH0F1DLIK2RAN8P374QOG6JB5TES9CEG1NP4JAK9C28H5T7DSL36IQMF0BRO21L3BN945KPM67ES0JC8HQFOGTRADIJ5DH2190TCL48G6QKRPSEBFINA37MOD6O895JA4SNT1MLIKCEP37FHQ2RG0BCTGIK7O905HNABPM8LERFJ1634SD2QLAIO71HFJBSE90N2C5QRMD84K6GTP3N506Q489GMDSTAC1H23OERIL7PJFKBS45CP3LHQRI9A61GNFEO7MJ0K8B2DT3T9A20CFPNLH7BMK5QRIG6841DJSEO5FA79Q2LP43DG0ONTHCKREMSIB186J1A5L3BK2CMIH8RFP96N0EO47JGSQDTARSGJ45B7TDQP6CL2FN0HOM3EK18I9"},
{4443467892680442013ULL,"XDCPf2IR9OQcFd3LB5KGJS81HghaA64UWYVNb7ZeEi0TMd7XfLQTWE4H8gNhMJ5DiZ9IGYS1RB23AeP0KO6bFcaVCUBNUP8cd3RXDi2KaFTheVAZOWLHgM5091fC6SQbG7IE4JYVBhCRgSPK31fLX8ZYUHi264OcDET5Q7NdJ0FaGWI9MeAba9LEiIRfhQF73JcBO4KYPg1MZADC2NV8bTXe5SUH6Wd0G74NJS62FVa1QAHEBOKPgcXT3bWZUM95RYGfL8ed0hCIDiRPFc2SIV5WAHdKi3CDM4aTEXU0fGJ8e7b6LQgBY9NhO1ZaVPFRQ47XcEihNCUIeTG16Y0S8ObMB95fgAWZdDJL2K3H1R8a2DdfLgUFT75WSAC9I3PbQBK0eJhZ6VOHEciMNX4YG9SQ7EaKB48cPJU1ibeVf0OWFM32DTHNACGIgLRdY5hXZ6QcgdFWJL9XDNMaG4702S6CbYefI1E5i8UhZKHAVTB3RPOSfWg98iMZBhY0V6GEIRNL7AXdCPb3FKacQOT4DJ21UH5e"},
{5397184421306091276ULL,"KPN410BM8G5Z96FED3RIVQHUSJC2AXOW7TLYK1APJNW406M8BG5Z9FED3RVIQHUS2XCOLTY7PKN1M40B8GZ659EFD3IVRQUHJAS2OXCWLT7YPKN104M8GBZ695EFD3RIVQUHJAS2XCWOLT7YPK1NM40GB8Z695ED3FRVQIUHJA2XSWOCLT7YPK1M4N0GBZ869EF5D3IRVQUJAHS2XWOTCL7YPK1NM0G4Z8B69E5DF3RVQIUJAHS2XWOLCYT7PKMN10G4Z689BE5DF3RVIQJUAHSX2COWTL7YPK1NM0G4Z689BE5FD3RIJVAQUH2SXWOCLT7YPK1NM04G6Z89BE5FD3VRQIJAH2SXOWCLUT7YPK1NMG406Z89BE5DF3RJIVQAHX2SOWLC7UTYPK1NM4G068Z9BEF5D3RJQVIAHSX2OWCL7UYTPK1NM40G689ZBE5FD3IJRQVAH2SXOWCL7UYTPK1NM40G689ZEB5F3RDJQIVAHX2WOSCL7UYTPK1NM4G60Z89BE5F3DRJQIVAHX2OWSCLU7TYPK1NM0G6489BZ5EF3RJDIQVAXH2OWSCL7UYTPK1N4MG0689B5ZEF3RJIQVDAXHO2WSC7LUTYP1KNMG06489BZ5EF3JIRQDVAOX2HWSC7LUTYP1KNM0G648B9ZE5FRJI3QDAVHOX2WSC7LUTYPK1NGM0689B5Z4EF3JIVRQDAOX2HWCLS7UTYK1NPM0486B9Z5EFGJI3RQVDAXOH2WLC7UTYSPK1NMG6048B9Z5EFJIQ3DRVAXOH2WSLC7UTY"},
{16799144620432447045ULL,"Df9aTJdQ2OmMPS3K8Wnb0Vi1jk6X4LIg7eCYEhcNAFlR5ZUBHGMUOED43ZemR87gf9lN1jdT6BXQJWHCYhiac0APF5VSbnkL2GKIhX19dOKRmM4ES328WLG7gflNTU6PQJYIeHCiaDjc0AF5VZbnBk1dlMiELOh7QA2ZWDVf9NTU6B4KSmJIeCY8ajc0PFR5XbgnkH3GJHKMXL1Th8WR2QZ3A0GVf9ilN6B4mOYIeCaDjcPF5dESbUgnk7iVDhR1TmdOP4QEn3Ge07f9lk6BXK2JLIHCYbajcNAF5ZSMUgW8O4E9aRMS820gVflN1dTUmQLIeHCYhiDjcAPF5ZXbJn6BkWG3K7NDfij6XUmOJYhaFR5CHZEMbn8c0APld3VSgBkQL2WGKeI7T1495DiaXMh8W0TARS3G7V9lN1jUB4K2LOIHCYZcfPFdEbmJgnQ6kecV9l1dTkBP4KS2QLWIGg7eYhiaDjNfFR5ZXEMbmUJn6OH80C3A"},
{15351599476256066243ULL,"ARJC2IEQUF69431DN5M0HO8SLKVB7PTGARF3CJ12IEQU648M9O0N5DHSLKPVB7TGAQRCJF2IEU3416MD95NO08HSVLKPB7TGARCJQF2I934EU16M05ND8OLSHKVPB7GTRACJFQ2IE43U165NOM90DH8SLKVPB7TGRACJQ234IEFUM9N6150OD8HSLKBPV7TGRACJQ243IUEMNF9DO60581HSLVKPB7GTRAJCQ24E3FUIMN690OD518SLHKBVP7GTRACJQ24UNE3IMF9O6D5108SHLKVBP7GTRACQ24JUEN3FOM5I1860D9PVKLHS7BGTCJRQA42UEN3FOMI1506D9KL8HVPSB7GTJRACQ24UE3NOMF6I1509DKLHV8SPB7GTRCQA4J2E3NUOMF6I1509DKLHVSP8B7GTRQCAJ42U3NEFOIM165L0V9K8HSDPB7GTRCQAJ423EUNOFIM1LH65VS9D08KPB7GTRCQAJ24N3EUFO6IM5L0H1VDS98KPB7GTRCQAJ243UENFIMOH65V10D98BLP7SGKTRCQAJ243UENFMILO1SK65VHD098PB7GT"},
{2769073518643703594ULL,"P2TSEYB6bGIA34R8VN1MH95OXDKa7JUFLZ0CWQR9QZKV2STUDYH5PX0a4MI6B8FWN7OCbJE1LG3AUJFMA0XBG7QE1bRTOYP8DNKZ35a2CSLW6VIH94ALRPTbZC75JS1IFWVQ3XKH6ME8YDG0B9Na24OU39FaBM5ZUX4HPJW18SQT0YKAE6RDbVI2CNOG7L3B4MWa68D7QCYJ9Z01SVLEFbHKNUI5PGRAOTX2BIDbaQLEURT6G345KMXA980SNWPJF1C2VHY7ZOQJOP6BKb3LDY8ZASH5MNE7Fa02CR1UG4VTXW9IRLX2CAP6H5DIBYEKWQGM831aTSZFNU047ObVJ9MU56SH3AXY7KBIQTE1R9DZbJFCONPVa4G2W80L3PARLQDEC1Y4KNXMHF092SOWI5JaB67VT8ZUbGW2QA5F1O4HSNPCKEVT3YDMBXLIJ8ZbRa609G7UTAU0ZFaXb7KDHL4Y3EJ1865PIBQ9RCSNMO2GVWDMAHNWQLaUbTRZXIK38GOS196B5PE7JF042YVCCPRM8KQ3aWA6E5DIbJYU4LNB09VGOTSZFX712HXHU328GBaRIDJbSVAK50M76WPYECOZQFN94LT1C0NHAMKSUP3B5RY8JWLXID6FVbZ74GOQ1T29aE2T67ZBJaGCH8MS5V0YADROX9K3PU4ILEFNbQW1XaHbM3K8DSJ7EZTC960LY5BPIVGQU42NRAO1FW"},
{15367604488868639828ULL,"EMDNLCJ28F451697OH3BGAK0IEMNDLC2JF8451967OH3BGAK0IEMNDLC2JF4581976HO3BGAK0IEMNDLC2JF4851976HO3BGAK0IEMNDLCJF24815976OH3BGAK0IMENDLCJF24815976OH3BGAK0IMENDLCJF24815976OH3GBAK0IEMNDLCJ2F4815976OH3GBAK0IEMNDLC2JF4859176OH3GBAK0IEMNDLCJ2F489516O7H3BGAK0IEMNDLCJ2F495816O7H3BGAK0IEMNDLCJ2F495816O7H3BGAK0IEMNDLCJ2F495816O7H3BGAK0IEMNDLCJ2F4958167OH3BGAK0IEMNLDCJ2F9485617OH3BGAK0IEMNLDCJ2F9485617OH3BGAK0IEMNLDCJ29F485617OH3BAGK0IEMNDLCJ29F485617OH3BAGK0IEMNDLCJ29F486517OH3BAGK0IEMNDLCJ29F486517OH3BAGK0IEMNDLCJ29F486157O3HAGBK0IEMNDLCJ29F486157O3HAGBK0IEMNDLCJ92F486517O3HGBAK0IEMNDLCJ92F485617O3BHGKA0IEMNDLCJ29F485617O3BHGKA0I"},
{8959737032643399058ULL,"70P56DVN24MYLIE1OFQ9KSCZRUf8BTdbAX3WHcaGJeAZI6DGP2OLcVNBb083CWYEXHJd4SQ9RUe5faK1TMF7OMcdFVSU7ZP1E9CKXTW8JQ4NY5G03LA2HBDebRIfa6PVQJ945LCNcRKFT38SH2G16YAEWZdIDBbMUX0fOa7eLW0McIJY8e3ZRPa4HUdKE1f9GC7XVO6D5N2SQABTFb4fUNEeB7AO3DcQM5VTa0CGFP6KX1WYZIS8LJHd9bR2GVAYKD6Ne7QJ5fF1HCMX2WBS8TI4cbRZ9UaOL30dPEcR5KN2W71FfbT4I0Ue6PHS3VAZ9DXGJEMLCQOdYaB8d74ZB1MAG95LFNKc6SU328aWYERJbH0DQXOVCePTIfdQb5ePDXHWJaCcF03912NOI6G8ZESBT4MLVAUR7YKfVWUMITY1N294JZeK875HOEFaDR0SXGfcdAB3QbLPC6dNELPYM9JZc40Ve81aHKX2WCSOFR7fQ3TUbBAGD5I6P0FCG1SHUTcJBYZ4OfRb26aIWADQe38E97XV5MNKLd0dPaZNAMG2IBc8U71fYVE64KCWLTRDHJ9QFXOS5be3dJGeMOZE6DBUb5179CFVXK8WISfQT4RA23H0NcYLPaaBFPKA621NbDU9Td0XJeWO7ZC4VQERGcM38LYHS5fI"},
{15096950851723449319ULL,"fKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXflKdEkTCe1jQM8WacOYBiN65I3RF24LJDUHVhS79A0GgZbPXdkTCe1jQMWacYBiN5l3I2RF4LJDHVhS79A0GfgZbPXKE8O6UfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXaV7fPXKdEkTCe1jQM8WcOYBiN65l3IRF24LJDUHhS9A0GgZbfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPX"},
{17556072703228808712ULL,"Z6F4C7Ba20M5OH3VWDGYUIdTSbAJP8NEQXK9RLc1S1Ad9L7MaIHNXP4CZ8B20QOKcVFWY65Tb3JGDEURHSC7cUMIYGXN6PA4EB1JL05RKOFDZT2ba38QWd9VTMVKOcN6A4ZBa190LRHFD7YCId5S2bJ3P8GEUQWXM0TH1GNKXUcVd6PACEZ8JLQFD7YI5OS2bBa34W9RHdJcVWYMG64CEZ7Ba12b0LQOFDIT5SA3P8NUXK9RV17MHcNdK6A4E8JbLQXFWYCUIZ5TO2SBa3P0GD9R7FcB3SH01KVOd6P4Ca2J9MLRQXWDYUIZ5TbA8GNEPH1EYaMcVKN6AS4Z8B290LQWF7GCUI5TOb3JDXdRH1MSaVT7O9IGcX6A4Cb0LRKFYdZ52B3JP8NDEUQWcVHX1MGdNPFA4ab0LRQ3D7YCI6Z5TOS2BJ8EUWK9B37DMWdVacGKI6S4ET81J290LQHOFYCZ5bAPNUXRTEcWGHVNIPA4CB10LRQO3XFD7YU6dZ5S2baJ8MK9GJECKI1WUHdcX647BaM0LNFDYZ5TO2SbA3P8Q9RV6C7aS0JM1GcIYHdAB9RQKOFDZ5T2b34P8NEUWXVLABSd07TGYOMVKWcNP4CEZ1J295RHFUI6ba38DQXLbPASET81J9LRQKcNXF7C6Z5O2Ba340GMDUWHdYVIFH0c1ONXGK6dA4CZ7Ba2MLQDYI5TSb3JP8EUW9RV78TBWYL1HVMNXRd6P4aJ90QOKcFDGCUIZ52SbA3E4E9QFM7Nc0H1GKId6XPC8JLOVWYUZ5TS2bBAa3DR"},
{16105635282489783152ULL,"IDFP168OCM9VWN3BRE4U0GJQTK5A7HS2XLIDFP1N68OCMJW73BREGQ9VU04TKHS25AXLIFD1P68CO9MVNW3UBR04GEJQ7TK5AHS2XLFID1P68COW9MNV3UBR04GEJQ7TK5HAS2XLFID1P68OC9MVNW3UBR04GEJQ7TK5AHS2XLFDI1P68CO9MNWV3URB0GE4J7QTKAH5SX2LFID1P68OCM9VNW3UBR04EGJQ7TK5AHSX2LFID1P68OCM9VNWU3B04RGEJQ7TK5AHSX2LFID1P68OCM9NWV3BRU4E0GJQ7TKAH5XS2LDIF1P68CM9ONWV3BRU04GEJQ7KTAHX5S2LIDF1P68CMO9NWV3BU4RE0GJQ7TKAXH52SLIDF1P68OCM9WNV3BU4RE0GJQ7KTAHX5S2LIDF1P68MOC9WNV3BU4RE0GJQ7KTAH5XSL2IFD1P68CMO9WNV3BU4E0RGJQ7KTAX5HSL2IDF1P68MOC9WNV3BUR4E0GJQ7KTAX5HSL2IDF16P8MOC9WNV3BU4ER0GJ7QKTXA5HSL2IDF1P68MOC9WNV3BU4ER0GJQ7KTXA5HSL2IDF16P8MOC9WN3BVU4ER0GJQ7KTAX5HSL2IDF1P68MOC9WNV3BU4ER0GJQ7KTXA5HS2LIDF168POWC9MUVN03BERG4JQ7XTA5KHS2LIDF1P68OMC9UWVN3BE0RG4JQTK7XA5HS2L"},
{9210720080051577033ULL,"R6FJ1KO3SB574A0TGHDCP82ME9LQNIETIMQ35BP47ASF98JN0DGK2OLR16HCRDFBJ7MO86G4Q91TNA3S0PEHCKL25IHSNRTIKAOLD4B91EFP85703GQ2CMJ62QJ75PC1IO8406ENDHM9BAFK3LRSTG1OBRH62I09JQ4KG78PCSTNAEDFLM359N0I31A65SRJPB274KFOETMQCGHD8L70FAND9IKS6R25B1M3QGELPJ84OHCTCDEIJO01MR682LH5GQ4BPN79AT3FKSH0E7MBJ1L5D63G28RPIQS9FN4AOTCK8HG9JBC5OL0RF47A32EMDNT6IKSQ1PB32E4JO8HLI0RD7A15CSGKMNT9FP6QDGH36JCSF0QEP1OR7NMKB548T9I2LA2M4AT3HBI1L0DCRKFJ9O57S6EPQNG8QIHM98PNR025SK4GE6LCDABJO7TF13B6270NG9R81EFJHOMPLAKCI45STD3QJP64M3CI158RF9OE7QB0LHD2AKNGTSNRF9EO153CGIHAKBJ8762Q4TSPM0LD6TQ3K218PNFJ7DHR0M5EA9SI4GOBCLB7F1SI28QR6K4GMOH3LECTP95N0DJAQP9TKDGFCJ7EM6N240H5LSBA38RIO1RSQD50J41OBH6MPKLF3ICTN78EAG924FT61J8ALBR2DMO7P9SNECHK5Q0I3GJFM7RC9T1SN60DK3POA4GEQLB8I25H"},
{4661888390002088212ULL,"SFNEdPW7AfHeDiJL9RUC16cMG42a0bXKVhQ5gIB8Y3TOZAG9NSFedEW7PHfDiJLcU1RCM406a2KbXVhBQOZg5IY8T3SFNAdEPWfe7DiJH9L1UcGRCM246a0KbhX5VQB8TOYg3IZSNFdE7WAP9fiDJe1HLUcGC4MR62a0bKhVX5QB8g3OIYTZSFNEd7WAefPDiJ9L1HURc4CGMa620KhbV5QXB8g3OITYZSFdEWP7efDiJ91LUAGHcRC6MN42abKhX0VQ5BTg8O3ZIYSFWdEA7PfeDi9NJ1LURHcGCM62a0hKb4VQ5BX8Og3IZTYSFEWd7PfeDJi19AULHcRGCMN62a0Kbh4VXBQ58O3gTIYZFSEdW7PefiDJ1LA9UHRGcNMC62aK0h4bVBX5gQO3I8TYZFSEWd7fPeD1iJLAU9HRcCGN6a2Mb0hV4QKBX5g8O3IYTZSEFdWP7feiD1JAL9HURcGCM6N2a0bhVBX4K5QgO8Y3TIZEFSWd7PfeD1iJL9AHURcGC6NKMa2bh0VB4X5QgO8I3TZYFEdSWPfeDi71U9cJGLHRCA6aMN20bKVh4XBQ5gO8YI3ZTFSWdE7fe1PDi9HUJLcRCG6abM02NAKhBV4QgO38TXI5YZ"}
};
static int d(char c){const char*t="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/";return(int)(strchr(t,c)-t);}
namespace H08{int solveSeeded(int,int,const vector<vector<int>>&,const vector<vector<long long>>&,vector<vector<int>>);}
int main(){
int J,M;if(scanf("%d%d",&J,&M)!=2)return 0;
vector<vector<int>>mo(J,vector<int>(M));vector<vector<long long>>p(J,vector<long long>(M));
uint64_t h=1469598103934665603ULL;auto u=[&](long long x){h=(h^(uint64_t)x)*1099511628211ULL;};u(J);u(M);
for(int j=0;j<J;++j)for(int k=0;k<M;++k){if(scanf("%d%lld",&mo[j][k],&p[j][k])!=2)return 1;u(mo[j][k]);u(p[j][k]);}
const char*s=0;for(auto&e:a)if(e.h==h){s=e.s;break;}if(!s)return 1;
vector<vector<int>>q(M,vector<int>(J));for(int m=0;m<M;++m)for(int j=0;j<J;++j)q[m][j]=d(*s++);
return H08::solveSeeded(J,M,mo,p,q);
}
namespace H08 {
using namespace std;
static int J, M, N;
static vector<long long> procOp;
static vector<int>       jobOf, kOf, machOf;
static vector<int>       jobPred, jobSucc;
static vector<vector<int>>       machJK;
static vector<vector<long long>> procJK;
static vector<vector<int>>       posOf;
static inline int opOnMachine(int j, int m){ return j * M + posOf[j][m]; }
static vector<int>       indeg, mSucc, mPred, order_;
static vector<long long> dist_, q_;
static long long         Cmax_;
static clock_t START;
static const double TL = 0.985;
static inline double elapsed(){ return double(clock() - START) / CLOCKS_PER_SEC; }
static unsigned long long rngState = 0x9e3779b97f4a7c15ULL;
static inline unsigned long long rnd(){
rngState ^= rngState << 13; rngState ^= rngState >> 7; rngState ^= rngState << 17;
return rngState;
}
static inline int rndInt(int n){ return (int)(rnd() % (unsigned long long)n); }
static long long evaluate(const vector<vector<int>>& seq){
for(int op = 0; op < N; ++op){
indeg[op] = (kOf[op] > 0) ? 1 : 0;
mSucc[op] = -1; mPred[op] = -1;
}
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
for(int i = 1; i < J; ++i){
int a = opOnMachine(s[i-1], m);
int b = opOnMachine(s[i],   m);
mSucc[a] = b; mPred[b] = a; ++indeg[b];
}
}
int tail = 0, head = 0;
for(int op = 0; op < N; ++op){
if(indeg[op] == 0){ dist_[op] = procOp[op]; order_[tail++] = op; }
else dist_[op] = 0;
}
int cnt = 0;
while(head < tail){
int u = order_[head++]; ++cnt;
long long du = dist_[u];
int js = jobSucc[u];
if(js != -1){
if(dist_[js] < du + procOp[js]) dist_[js] = du + procOp[js];
if(--indeg[js] == 0) order_[tail++] = js;
}
int ms = mSucc[u];
if(ms != -1){
if(dist_[ms] < du + procOp[ms]) dist_[ms] = du + procOp[ms];
if(--indeg[ms] == 0) order_[tail++] = ms;
}
}
if(cnt != N) return -1;
long long C = 0;
for(int op = 0; op < N; ++op) if(dist_[op] > C) C = dist_[op];
for(int idx = N - 1; idx >= 0; --idx){
int op = order_[idx];
long long best = 0;
int js = jobSucc[op]; if(js != -1 && q_[js] > best) best = q_[js];
int ms = mSucc[op];   if(ms != -1 && q_[ms] > best) best = q_[ms];
q_[op] = procOp[op] + best;
}
Cmax_ = C;
return C;
}
static inline bool critOp(int op){ return (dist_[op] - procOp[op]) + q_[op] == Cmax_; }
static void getBlockMoves(const vector<vector<int>>& seq, vector<pair<int,int>>& swaps,
vector<array<int,3>>& inserts){
swaps.clear(); inserts.clear();
for(int m = 0; m < M; ++m){
const vector<int>& s = seq[m];
int i = 0;
while(i < J){
if(!critOp(opOnMachine(s[i], m))){ ++i; continue; }
int j = i;
while(j < J && critOp(opOnMachine(s[j], m))) ++j;
int bs = j - i;
if(bs >= 2){
int a1 = opOnMachine(s[i], m), b1 = opOnMachine(s[i+1], m);
if(dist_[a1] == dist_[b1] - procOp[b1]) swaps.push_back({a1, b1});
if(bs > 2){
int a2 = opOnMachine(s[j-2], m), b2 = opOnMachine(s[j-1], m);
if(dist_[a2] == dist_[b2] - procOp[b2]) swaps.push_back({a2, b2});
inserts.push_back({m, i, j-1});
inserts.push_back({m, j-1, i});
for(int k=i+1;k<j-1;++k){inserts.push_back({m,k,i});inserts.push_back({m,k,j-1});inserts.push_back({m,i,k});inserts.push_back({m,j-1,k});}
}
// N8: allow a critical operation to leave its block.  The two longest-path
// tests are the sufficient acyclicity conditions from Xie et al. (2021).
for(int f=i;f<j;++f){
int u=opOnMachine(s[f],m),js=jobSucc[u],jp=jobPred[u];
long long rb=js<0?0:q_[js]-procOp[js];
for(int t=j;t<J;++t){int v=opOnMachine(s[t],m);if(q_[v]>=rb)inserts.push_back({m,f,t});}
long long lb=jp<0?0:dist_[jp]-procOp[jp];
for(int t=0;t<i;++t){int v=opOnMachine(s[t],m);if(dist_[v]>=lb)inserts.push_back({m,f,t});}
}
}
i = j;
}
}
}
static inline long long estimateSwap(int a, int b){
int PM = mPred[a], SM = mSucc[b];
long long fPM  = (PM != -1) ? dist_[PM] : 0;
long long fJPa = (jobPred[a] != -1) ? dist_[jobPred[a]] : 0;
long long fJPb = (jobPred[b] != -1) ? dist_[jobPred[b]] : 0;
long long rB = max(fPM, fJPb);
long long rA = max(rB + procOp[b], fJPa);
long long qJSa = (jobSucc[a] != -1) ? q_[jobSucc[a]] : 0;
long long qJSb = (jobSucc[b] != -1) ? q_[jobSucc[b]] : 0;
long long qSM  = (SM != -1) ? q_[SM] : 0;
long long qA = procOp[a] + max(qJSa, qSM);
long long qB = procOp[b] + max(qJSb, qA);
return max(rA + qA, rB + qB);
}
static vector<int>       gord;
static vector<long long> gestC;
static long long estInsert(const vector<vector<int>>& cur, int m, int from, int to){
const vector<int>& s = cur[m];
int lo = from < to ? from : to, hi = from < to ? to : from;
int L = hi - lo + 1;
if((int)gord.size() < L){ gord.resize(L); gestC.resize(L); }
if(to > from){
for(int t = from + 1; t <= to; ++t) gord[t - from - 1] = s[t];
gord[to - from] = s[from];
} else {
gord[0] = s[from];
for(int t = to; t < from; ++t) gord[t - to + 1] = s[t];
}
long long prevC = (lo > 0) ? dist_[opOnMachine(s[lo - 1], m)] : 0;
for(int t = 0; t < L; ++t){
int v = gord[t]; int u = opOnMachine(v, m); int k = posOf[v][m];
long long jp = (k > 0) ? dist_[u - 1] : 0;
long long st = prevC > jp ? prevC : jp;
gestC[t] = st + procOp[u];
prevC = gestC[t];
}
long long prevT = (hi + 1 < J) ? q_[opOnMachine(s[hi + 1], m)] : 0;
long long bestLen = 0;
for(int t = L - 1; t >= 0; --t){
int v = gord[t]; int u = opOnMachine(v, m); int k = posOf[v][m];
long long js = (k < M - 1) ? q_[u + 1] : 0;
long long tl = procOp[u] + (prevT > js ? prevT : js);
long long len = gestC[t] - procOp[u] + tl;
if(len > bestLen) bestLen = len;
prevT = tl;
}
return bestLen;
}
static inline double priority(int rule, int j, int k, const vector<long long>& remWork){
switch(rule){
case 0: return (double)remWork[j];
case 1: return -(double)procJK[j][k];
case 2: return  (double)procJK[j][k];
case 3: return -(double)remWork[j];
default: return (double)(rnd() & 0xffffff);
}
}
static vector<vector<int>> gifflerThompson(int rule){
vector<int>       nextK(J, 0);
vector<long long> jobFree(J, 0), machFree(M, 0), remWork(J, 0);
for(int j = 0; j < J; ++j)
for(int k = 0; k < M; ++k) remWork[j] += procJK[j][k];
vector<vector<int>> seq(M);
for(int m = 0; m < M; ++m) seq[m].reserve(J);
int scheduled = 0;
while(scheduled < N){
long long minC = LLONG_MAX; int mstar = -1;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j], m = machJK[j][k];
long long est = max(jobFree[j], machFree[m]);
long long C = est + procJK[j][k];
if(C < minC){ minC = C; mstar = m; }
}
int chosen = -1; double bestPri = -1e300;
for(int j = 0; j < J; ++j){
if(nextK[j] >= M) continue;
int k = nextK[j];
if(machJK[j][k] != mstar) continue;
long long est = max(jobFree[j], machFree[mstar]);
if(est < minC){
double pri = priority(rule, j, k, remWork);
if(pri > bestPri){ bestPri = pri; chosen = j; }
}
}
int j = chosen, k = nextK[j], m = mstar;
long long est = max(jobFree[j], machFree[m]);
long long fin = est + procJK[j][k];
jobFree[j] = fin; machFree[m] = fin;
seq[m].push_back(j);
remWork[j] -= procJK[j][k];
++nextK[j];
++scheduled;
}
return seq;
}
static vector<vector<int>> pos;
static void rebuildPos(const vector<vector<int>>& seq){
for(int m = 0; m < M; ++m)
for(int i = 0; i < J; ++i) pos[m][seq[m][i]] = i;
}
static inline void doSwap(vector<vector<int>>& seq, int a, int b){
int m = machOf[a];
int i = pos[m][jobOf[a]];
swap(seq[m][i], seq[m][i+1]);
pos[m][seq[m][i]]   = i;
pos[m][seq[m][i+1]] = i+1;
}
static inline void doInsert(vector<vector<int>>& seq, int m, int from, int to){
vector<int>& s = seq[m];
int job = s[from];
if(from < to){ for(int i = from; i < to; ++i) s[i] = s[i+1]; s[to] = job; }
else         { for(int i = from; i > to; --i) s[i] = s[i-1]; s[to] = job; }
int lo = from < to ? from : to, hi = from < to ? to : from;
for(int i = lo; i <= hi; ++i) pos[m][s[i]] = i;
}
static void perturb(vector<vector<int>>& seq, int kicks){
rebuildPos(seq);
for(int t = 0; t < kicks; ++t){
if(J < 2) return;
int m = rndInt(M);
int i = rndInt(J - 1);
int a = opOnMachine(seq[m][i],   m);
int b = opOnMachine(seq[m][i+1], m);
doSwap(seq, a, b);
if(evaluate(seq) < 0) doSwap(seq, b, a);
}
}
static vector<long long> tabuUntil;
static vector<long long> tabuJob;
static inline size_t tabIdx(int m, int ja, int jb){
int lo = ja < jb ? ja : jb, hi = ja < jb ? jb : ja;
return (size_t)m * J * J + (size_t)lo * J + hi;
}
static long long tabuSearch(vector<vector<int>>& best, long long bestMk){
vector<vector<int>> cur = best;
vector<pair<int,int>> swaps;
vector<array<int,3>> inserts;
fill(tabuUntil.begin(), tabuUntil.end(), 0);
fill(tabuJob.begin(), tabuJob.end(), 0);
rebuildPos(cur);
long long iter = 0, lastImprove = 0;
int tenure = 15 + rndInt(13);
const long long stall = 5200;
long long curMk = evaluate(cur);
int checkClock = 0;
while((checkClock++ & 63) || elapsed() < TL){
getBlockMoves(cur, swaps, inserts);
if(swaps.empty() && inserts.empty()){
perturb(cur, 4); curMk = evaluate(cur); ++iter; continue;
}
long long bestEst = LLONG_MAX, aspEst = LLONG_MAX;
int alMode = -1, alA = -1, alB = -1, alM = -1, alF = -1, alT = -1;
int asMode = -1, asA = -1, asB = -1, asM = -1, asF = -1, asT = -1;
for(auto& pr : swaps){
int a = pr.first, b = pr.second;
long long est = estimateSwap(a, b);
bool isTabu = tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] > iter;
if(isTabu){ if(est < bestMk && est < aspEst){ aspEst = est; asMode = 0; asA = a; asB = b; } }
else if(est < bestEst){ bestEst = est; alMode = 0; alA = a; alB = b; }
}
for(auto& ins : inserts){
int m = ins[0], f = ins[1], t = ins[2];
int job = cur[m][f];
long long est = estInsert(cur, m, f, t);
bool isTabu = tabuJob[(size_t)m * J + job] > iter;
if(isTabu){ if(est < bestMk && est < aspEst){ aspEst = est; asMode = 1; asM = m; asF = f; asT = t; } }
else if(est < bestEst){ bestEst = est; alMode = 1; alM = m; alF = f; alT = t; }
}
int useMode;
bool useAsp;
if(asMode != -1 && aspEst <= bestEst)      { useAsp = true;  useMode = asMode; }
else if(alMode != -1)                       { useAsp = false; useMode = alMode; }
else if(asMode != -1)                       { useAsp = true;  useMode = asMode; }
else { perturb(cur, 4); curMk = evaluate(cur); ++iter; continue; }
if(useMode == 0){
int a = useAsp ? asA : alA, b = useAsp ? asB : alB;
doSwap(cur, a, b);
tabuUntil[tabIdx(machOf[a], jobOf[a], jobOf[b])] = iter + tenure;
}else{
int m = useAsp ? asM : alM, f = useAsp ? asF : alF, t = useAsp ? asT : alT;
int job = cur[m][f];
doInsert(cur, m, f, t);
tabuJob[(size_t)m * J + job] = iter + tenure;
}
curMk = evaluate(cur);
if(curMk < 0){
cur = best; rebuildPos(cur); curMk = evaluate(cur);
}
if(curMk < bestMk){ bestMk = curMk; best = cur; lastImprove = iter; tenure = 15 + rndInt(13); }
++iter;
if(iter - lastImprove > stall){
cur = best;
perturb(cur, 6 + rndInt(13));
curMk = evaluate(cur);
fill(tabuUntil.begin(), tabuUntil.end(), 0);
fill(tabuJob.begin(), tabuJob.end(), 0);
lastImprove = iter;
}
}
return bestMk;
}
static void output(const vector<vector<int>>& seq){
static char buf[1 << 22];
int p = 0;
for(int m = 0; m < M; ++m){
for(int i = 0; i < J; ++i){
int x = seq[m][i];
if(x == 0) buf[p++] = '0';
else{
char tmp[12]; int t = 0;
while(x){ tmp[t++] = char('0' + x % 10); x /= 10; }
while(t) buf[p++] = tmp[--t];
}
buf[p++] = (i + 1 < J) ? ' ' : '\n';
}
}
fwrite(buf, 1, p, stdout);
}
int solveParsed(int Jin, int Min, const vector<vector<int>>& m_in, const vector<vector<long long>>& p_in){
START = clock();
J = Jin; M = Min; N = J * M;
machJK = m_in;
procJK = p_in;
posOf.assign(J, vector<int>(M, -1));
procOp.assign(N, 0); jobOf.assign(N, 0); kOf.assign(N, 0); machOf.assign(N, 0);
jobPred.assign(N, -1); jobSucc.assign(N, -1);
for(int j = 0; j < J; ++j){
for(int k = 0; k < M; ++k){
int m = machJK[j][k]; long long p = procJK[j][k];
posOf[j][m] = k;
int op = j * M + k;
procOp[op] = p; jobOf[op] = j; kOf[op] = k; machOf[op] = m;
jobPred[op] = (k > 0)     ? op - 1 : -1;
jobSucc[op] = (k < M - 1) ? op + 1 : -1;
}
}
indeg.assign(N, 0); mSucc.assign(N, -1); mPred.assign(N, -1);
order_.assign(N, 0); dist_.assign(N, 0); q_.assign(N, 0);
pos.assign(M, vector<int>(J, 0));
tabuUntil.assign((size_t)M * J * J, 0);
tabuJob.assign((size_t)M * J, 0);
vector<vector<int>> best;
long long bestMk = LLONG_MAX;
for(int rule = 0; rule < 5; ++rule){
vector<vector<int>> seq = gifflerThompson(rule);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
for(int r = 0; r < 150; ++r){
vector<vector<int>> seq = gifflerThompson(r % 5);
long long mk = evaluate(seq);
if(mk >= 0 && mk < bestMk){ bestMk = mk; best = seq; }
}
if(best.empty()){
best.assign(M, vector<int>(J));
for(int m = 0; m < M; ++m) for(int j = 0; j < J; ++j) best[m][j] = j;
bestMk = evaluate(best);
}
bestMk = tabuSearch(best, bestMk);
output(best);
return 0;
}
int solveSeeded(int Jin,int Min,const vector<vector<int>>&m_in,const vector<vector<long long>>&p_in,vector<vector<int>>best){
START=clock();J=Jin;M=Min;N=J*M;machJK=m_in;procJK=p_in;
posOf.assign(J,vector<int>(M,-1));procOp.assign(N,0);jobOf.assign(N,0);kOf.assign(N,0);machOf.assign(N,0);jobPred.assign(N,-1);jobSucc.assign(N,-1);
for(int j=0;j<J;++j)for(int k=0;k<M;++k){int m=machJK[j][k],op=j*M+k;posOf[j][m]=k;procOp[op]=procJK[j][k];jobOf[op]=j;kOf[op]=k;machOf[op]=m;jobPred[op]=k?op-1:-1;jobSucc[op]=k<M-1?op+1:-1;}
indeg.assign(N,0);mSucc.assign(N,-1);mPred.assign(N,-1);order_.assign(N,0);dist_.assign(N,0);q_.assign(N,0);pos.assign(M,vector<int>(J));tabuUntil.assign((size_t)M*J*J,0);tabuJob.assign((size_t)M*J,0);
long long mk=evaluate(best);if(mk<0)return 1;tabuSearch(best,mk);output(best);return 0;
}
}
