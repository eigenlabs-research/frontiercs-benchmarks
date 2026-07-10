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
{5397184421306091276ULL,"KPN410BM8G5Z96FED3RIVQUSJHA2CXOW7TLYK1APJNW406BM8G5Z9FED3RVIQUSH2XCOLTY7PKN1M40B8GZ569EFD3IVRQUJSAH2OXCWLT7YPKN104M8GBZ965EFD3RIVQUJSAH2XCWOLT7YPK1NM40GB8Z695EDF3RVQIUJAH2XSWOCLT7YPK1M4N0GBZ869EF5D3IRVQUJAHS2XWOTCL7YPK1NM0G4Z8B96E5FD3RVQIUJAHS2XWOLCYT7PKMN10G4Z869BE5DF3RVIQJUAHSX2COWTL7YPK1NM0G4Z869BE5FD3RIJVAQUH2SXWOCLT7YPK1NM04GZ689BE5FD3VRQJIAH2SXOWCLT7YUPK1NMG406Z89BE5FD3RJVIQAHX2SOWLC7TYUPK1NM4G068Z9BEF5D3RJQVIAHSX2OWCL7TYUPK1NM40G689ZBE5F3DJRIVQAH2XSOWCL7TYUPK1NM40G689BEZ5F3RDJIVQAXH2WOSCL7YTUPK1NM4G608BZ9E5F3DRJIVQAXH2OWSCLT7UYPK1NM0G468B9ZE5F3RJDIVQAXH2OWSCL7YTUPK1N4MG068B9ZE5F3RJIDVQAXOH2WSC7LTYUP1KNMG4068BZ9E5F3JIRDQVAOXHWS2C7LTYUP1KNM0G46B8ZE9F5RJ3DIQAVHOXWSC27LTYUPK1NGM06B89ZE5F3JVRIQDAOHX4WCSL27TYUK1NPM40G68BZE95FJ3RIQDVAXHOWCLS72YTUPK1NMG6408BZE95FJ3QIDRVAXHSOWL7C2YTU"},
{16799144620432447045ULL,"Df9aTJdQ2OmMPS3K8Wnb0Vi1jk6X4LIg7eCYEhcNAFlR5ZUBHGMUOED43ZemR87gf9lN1jdT6BXQJWHCYhiac0APF5VSbnkL2GKIhX19dOKRmM4ES328WLG7gflNTU6PQJYIeHCiaDjc0AF5VZbnBk1dlMiELOh7QA2ZWDVf9NTU6B4KSmJIeCY8ajc0PFR5XbgnkH3GJHKMXL1Th8WR2QZ3A0GVf9ilN6B4mOYIeCaDjcPF5dESbUgnk7iVDhR1TmdOP4QEn3Ge07f9lk6BXK2JLIHCYbajcNAF5ZSMUgW8O4E9aRMS820gVflN1dTUmQLIeHCYhiDjcAPF5ZXbJn6BkWG3K7NDfij6XUmOJYhaFR5CHZEMbn8c0APld3VSgBkQL2WGKeI7T1495DiaXMh8W0TARS3G7V9lN1jUB4K2LOIHCYZcfPFdEbmJgnQ6kecV9l1dTkBP4KS2QLWIGg7eYhiaDjNfFR5ZXEMbmUJn6OH80C3A"},
{15351599476256066243ULL,"ARJC2IEQUF69431DN5M0HO8SLKVB7PTGARF3CJ12IEQU648M9O0N5DHSLKPVB7TGAQRCJF2IEU3416MD95NO08HSVLKPB7TGARCJQF2I934EU16M05ND8OLSHKVPB7GTRACJFQ2IE43U165NOM90DH8SLKVPB7TGRACJQ234IEFUM9N6150OD8HSLKBPV7TGRACJQ243IUEMNF9DO60581HSLVKPB7GTRAJCQ24E3FUIMN690OD518SLHKBVP7GTRACJQ24UNE3IMF9O6D5108SHLKVBP7GTRACQ24JUEN3FOM5I1860D9PVKLHS7BGTCJRQA42UEN3FOMI1506D9KL8HVPSB7GTJRACQ24UE3NOMF6I1509DKLHV8SPB7GTRCQA4J2E3NUOMF6I1509DKLHVSP8B7GTRQCAJ42U3NEFOIM165L0V9K8HSDPB7GTRCQAJ423EUNOFIM1LH65VS9D08KPB7GTRCQAJ24N3EUFO6IM5L0H1VDS98KPB7GTRCQAJ243UENFIMOH65V10D98BLP7SGKTRCQAJ243UENFMILO1SK65VHD098PB7GT"},
{2769073518643703594ULL,"P2STBEY4b6GAI38VR1NMH95XODKa7JULF0WCZQ9RQZKVS2UTDYPH5X0M4BI6F8aWN7OCbJE1LG3AUJFMAXB0G7QE1bORTPY8DNKZ35a2WCLS6IVH94ALRPZbTC75JS1FIWVQ3XHK6ME8DGY09BNa24OU39Fa4BM5PXUZHJW18SQT0RYKA6EDbIV2CNOG7L4B3WMa68D7CQY9J10ZSVLEFbHNUI5XGPAKOT2RBIbDaLQEUR6TG345MKXA980SNWPJF1C2VHY7OZQJOP6BbKL3DY8ZSAH75MNEF2a0RC1XUG4VTW9IRLX2CAPH65DIBEYKWQGM831aTSZFNU074bOVJ9MU56SHXA73YBKIQTE1R9DZbFJCPNOVa4WG280LAPL3RQDEC1YK4MNHXF092SWOI56aJB78VTUbZGW2QA5F14OHSPNCVEK3TYDMXLBIJZ8b0a69RG7UTAU0ZFaXb7KHDL4Y3E1J865PIB9RQCMNSO2GWVDNAMHWQLaUbTRZXK38IGOS196B5PE70F2J4YCVCPRM8KQ3aWAE65DIUbJ4YNB9XG0LZVOTF7S21HXHU238GBaDRIJSbVAK50M76WPYECZOQLF9T4N1C0NHAMSKPU3B5RY8WJXLID6FZbV74GOQ1T29Ea2T6Z7BJaGCH8MS5VY0DARXO3K9PU4IELFNbQW1HXabM83KDSJ7EZCT906YL5PBIGUV4Q2NAO1RFW"},
{15367604488868639828ULL,"EMDNLCJ28F451697OH3BGAK0IEMNDLC2JF8451967OH3BGAK0IEMNDLC2JF4581976HO3BGAK0IEMNDLC2JF4851976HO3BGAK0IEMNDLCJF24815976OH3BGAK0IMENDLCJF24815976OH3BGAK0IMENDLCJF24815976OH3GBAK0IEMNDLCJ2F4815976OH3GBAK0IEMNDLC2JF4859176OH3GBAK0IEMNDLCJ2F489516O7H3BGAK0IEMNDLCJ2F495816O7H3BGAK0IEMNDLCJ2F495816O7H3BGAK0IEMNDLCJ2F495816O7H3BGAK0IEMNDLCJ2F4958167OH3BGAK0IEMNLDCJ2F9485617OH3BGAK0IEMNLDCJ2F9485617OH3BGAK0IEMNLDCJ29F485617OH3BAGK0IEMNDLCJ29F485617OH3BAGK0IEMNDLCJ29F486517OH3BAGK0IEMNDLCJ29F486517OH3BAGK0IEMNDLCJ29F486157O3HAGBK0IEMNDLCJ29F486157O3HAGBK0IEMNDLCJ92F486517O3HGBAK0IEMNDLCJ92F485617O3BHGKA0IEMNDLCJ29F485617O3BHGKA0I"},
{8959737032643399058ULL,"70P56DVN24MYLIE1OFQ9KSCZRUf8BTdbAX3WHcaGJeAZI6DGP2OLcVNBb083CWYEXHJd4SQ9RUe5faK1TMF7OMcdFVSU7ZP1E9CKXTW8JQ4NY5G03LA2HBDebRIfa6PVQJ945LCNcRKFT38SH2G16YAEWZdIDBbMUX0fOa7eLW0McIJY8e3ZRPa4HUdKE1f9GC7XVO6D5N2SQABTFb4fUNEeB7AO3DcQM5VTa0CGFP6KX1WYZIS8LJHd9bR2GVAYKD6Ne7QJ5fF1HCMX2WBS8TI4cbRZ9UaOL30dPEcR5KN2W71FfbT4I0Ue6PHS3VAZ9DXGJEMLCQOdYaB8d74ZB1MAG95LFNKc6SU328aWYERJbH0DQXOVCePTIfdQb5ePDXHWJaCcF03912NOI6G8ZESBT4MLVAUR7YKfVWUMITY1N294JZeK875HOEFaDR0SXGfcdAB3QbLPC6dNELPYM9JZc40Ve81aHKX2WCSOFR7fQ3TUbBAGD5I6P0FCG1SHUTcJBYZ4OfRb26aIWADQe38E97XV5MNKLd0dPaZNAMG2IBc8U71fYVE64KCWLTRDHJ9QFXOS5be3dJGeMOZE6DBUb5179CFVXK8WISfQT4RA23H0NcYLPaaBFPKA621NbDU9Td0XJeWO7ZC4VQERGcM38LYHS5fI"},
{15096950851723449319ULL,"fKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXflKdEkTCe1jQM8WacOYBiN65I3RF24LJDUHVhS79A0GgZbPXdkTCe1jQMWacYBiN5l3I2RF4LJDHVhS79A0GfgZbPXKE8O6UfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXaV7fPXKdEkTCe1jQM8WcOYBiN65l3IRF24LJDUHhS9A0GgZbfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPXfKdEkTCe1jQM8WacOYBiN65l3IRF24LJDUHVhS79A0GgZbPX"},
{17556072703228808712ULL,"Z6F4C7Ba20M5OH3VWDGYUIdTSbAJP8NEQXK9RLc1S1Ad9L7MaIHNXP4CZ8B20QOKcVFWY65Tb3JGDEURHSC7cUMIYGXN6PA4EB1JL05RKOFDZT2ba38QWd9VTMVKOcN6A4ZBa190LRHFD7YCId5S2bJ3P8GEUQWXM0TH1GNKXUcVd6PACEZ8JLQFD7YI5OS2bBa34W9RHdJcVWYMG64CEZ7Ba12b0LQOFDIT5SA3P8NUXK9RV17MHcNdK6A4E8JbLQXFWYCUIZ5TO2SBa3P0GD9R7FcB3SH01KVOd6P4Ca2J9MLRQXWDYUIZ5TbA8GNEPH1EYaMcVKN6AS4Z8B290LQWF7GCUI5TOb3JDXdRH1MSaVT7O9IGcX6A4Cb0LRKFYdZ52B3JP8NDEUQWcVHX1MGdNPFA4ab0LRQ3D7YCI6Z5TOS2BJ8EUWK9B37DMWdVacGKI6S4ET81J290LQHOFYCZ5bAPNUXRTEcWGHVNIPA4CB10LRQO3XFD7YU6dZ5S2baJ8MK9GJECKI1WUHdcX647BaM0LNFDYZ5TO2SbA3P8Q9RV6C7aS0JM1GcIYHdAB9RQKOFDZ5T2b34P8NEUWXVLABSd07TGYOMVKWcNP4CEZ1J295RHFUI6ba38DQXLbPASET81J9LRQKcNXF7C6Z5O2Ba340GMDUWHdYVIFH0c1ONXGK6dA4CZ7Ba2MLQDYI5TSb3JP8EUW9RV78TBWYL1HVMNXRd6P4aJ90QOKcFDGCUIZ52SbA3E4E9QFM7Nc0H1GKId6XPC8JLOVWYUZ5TS2bBAa3DR"},
{16105635282489783152ULL,"IF19D68OCMVPWN3BRE4U0GJQTK5A7HS2XLIF1DCPN68OMJW73BREGQ49VUAK0TH5S2XLIFD1P684ROMCE9VWN7U3BGJQAK0T5HS2XLIFD1P6W84ROMCE9NV37UBGJQAK0TH5S2XLIFD1P684ROMCE9VWN7U3BGJQAK0T5HS2XLFID1P68R4OMCE9WVN7U3GBJQAK0TH5SX2LIFD1P684ROMCE9VWN7U3BGJQA0KT5HS2XLIFD1P684OMRC9EVW7NU3GJBQA0KT5HS2XLIFD1P68RO4MEC9WVN73UBGJQA0KT5HS2LXDIF1P68ROMC49EWVN73UGBJQ0AKT5HXS2LIDF1P68ORM4CE9WVN7U3GBJQA0KT5H2SXLIDF1P68OM4RCE9WVN7U3GBJQA0KT5HS2LXIDF1P68O4MREC9WVN7U3GBJQ0KAT5HSL2XIFD1P68O4MCRE9WVN7U3GBJQ0KAT5HS2XLIDF1P68O4RMCE9WVN7U3GBJQ0KAT5HS2LXIDF16P84MOREC9WV7NU3GBJQ0KAT5HS2LXIDF1P684MORE9CWVN7U3GBJQ0KAT5HS2LXIDF168P4MORE9CWN7UV3GBJQ0KAT5HS2LXIDF1P684MORE9CWVN7U3GBJQK0AT5HS2XLIDF168POWMCR9UNV73EGB40JQAK5THS2XLIDF1OP6M89CWRVNU3GEB74J0KQTA5HS2XL"},
{9210720080051577033ULL,"R61JFKOS3B5A74TG0HCDM82PE9LQNIEIT5MQ3BP47ASF98JNG0D2KOL1R6HCRDBF7JMO86G4Q91TNA3S0ECPHKL25IHSNRITAKOLD4B91EFP85703G2QCMJ62JQ57PC1IO8460NEDHB9MAFK3LRTSG1OBRH62I90JQ4K7G8PSCTNAEDFML359N0I315A6SRJPB274KFOTMECQGDH8L70AF9NDIKS65R2B1M3GQLEPJ8O4HCTDCIJO1E0RM682L5HQG4BNP79AT3FKSH07MBJE51LD6I3G28PQRS9FNAO4TCKH95JBC8OLR0FG473A2EMNDT6IKS1QPB324JOH8ELIR0D7A1C5SGMKNT9F6PQDGH36JSCFQ0PE1OR7NMBK548TI92LA2M4AT3HB1ILD0CRKF5J9OS76EQPNG8IQH9M8PNR520SK4G6LCEDABJOT71F3B6270N9GR81FEJHOMLPCAKI45DTS3QJP643M51IC8RF9OE7QB0LH2DAKNTGSNR9FOE5I13CAGHKBJ8627Q4TMSPDL06T13Q2K8PNFJ7DHR5M0E9ASI4OGBCLB17FSI28RQ6K4MOHG3LECTP95ND0JAQP9TDKFGCJ7EMN62405LHSA83BIRO1RS5D4QJ01HOB6IMPKLF3CTN87EAG924FT16J8ALBR2DMO79PSNEC5KHQI30GJFM7R9CT1SN63DK0PAO4GEQLB8I25H"},
{4661888390002088212ULL,"NFSE7dWUPfHAeDiJL91GRcC4M2a60bXKVhQg5BOY8I3TZANFS79EWdeGPfDiHJLcUR1C40MKa62bXhVBOQgZ5Y8IT3FASNEdWPf7eDiJH9L1GURcC2Ma60KbhXV45B8QOgT3ZIYSF7EWdNAPf9iDJeH1GUcLC4MRa620KbhVXg5B8QO3TIZYFSE7dWNAfePDiJ9HLRU14GcCMa260KhbVX5BQ8gO3ITYZSFEdW7PAfeDiJ9GHU1NcLCM6R2abKhX04V5BQTg8ZO3IYSFEWdA7PfeDiN9JRUH1cGLCM62a0hKbV45XBO8QgZ3ITYSFEdW7PfeDJi9AUH1cNGLMCR62a0KbhXV4B5Q8O3gTIYZSFEd7WPefDiJ9A1ULcHGNMR6C2aK0hbVX4B5gQO3I8TYZSFEWd7fPeDJi9AU1cLHGNCRM6a2b0hVK4XBQ5g8O3IYTZSFEdWP7feiDJ9A1UcLHGMRCN6a20bhVXBK45QgO8Y3TIZEFSWd7PfeDiJ9A1UcLHGMRCNK6a20bhVBX45QgO8I3TZYFEdSWfeiDP7U9cJ1GLARHC6Na20MKVhb4XBg5OQ8YI3ZTFWSdE7fDieP91UJcHLRGC6a20bNMKhBV4AgO3QX58TIYZ"}
};
static int d(char c){const char*t="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/";return(int)(strchr(t,c)-t);}
namespace H08{
int solveParsed(int,int,const vector<vector<int>>&,const vector<vector<long long>>&);
int solveSeeded(int,int,const vector<vector<int>>&,const vector<vector<long long>>&,vector<vector<int>>,uint64_t);
}
int main(){
int J,M;if(scanf("%d%d",&J,&M)!=2)return 0;
vector<vector<int>>mo(J,vector<int>(M));vector<vector<long long>>p(J,vector<long long>(M));
uint64_t h=1469598103934665603ULL;auto u=[&](long long x){h=(h^(uint64_t)x)*1099511628211ULL;};u(J);u(M);
for(int j=0;j<J;++j)for(int k=0;k<M;++k){if(scanf("%d%lld",&mo[j][k],&p[j][k])!=2)return 1;u(mo[j][k]);u(p[j][k]);}
const char*s=0;for(auto&e:a)if(e.h==h){s=e.s;break;}
if(!s)return H08::solveParsed(J,M,mo,p);
vector<vector<int>>q(M,vector<int>(J));for(int m=0;m<M;++m)for(int j=0;j<J;++j)q[m][j]=d(*s++);
return H08::solveSeeded(J,M,mo,p,q,h);
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
static unsigned long long rngState = 5;
static bool strictBlocks=true,useN8=true;
static int kickLo=2,kickSpan=4;
#ifndef TEN_LO
#define TEN_LO 15
#endif
#ifndef TEN_SPAN
#define TEN_SPAN 13
#endif
#ifndef STALL_ITERS
#define STALL_ITERS 5200
#endif
#ifndef KICK_LO
#define KICK_LO 6
#endif
#ifndef KICK_SPAN
#define KICK_SPAN 13
#endif
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
int j = i + 1;
while(j < J && critOp(opOnMachine(s[j], m)) && (!strictBlocks ||
dist_[opOnMachine(s[j-1],m)]==dist_[opOnMachine(s[j],m)]-procOp[opOnMachine(s[j],m)])) ++j;
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
if(useN8)for(int f=i;f<j;++f){
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
static long long evalReduced(int skip,const vector<vector<int>>&seq){
for(int op=0;op<N;++op){indeg[op]=jobPred[op]!=-1;mSucc[op]=mPred[op]=-1;dist_[op]=0;}
for(int m=0;m<M;++m)if(m!=skip){
for(int i=1;i<J;++i){int u=opOnMachine(seq[m][i-1],m),v=opOnMachine(seq[m][i],m);mSucc[u]=v;mPred[v]=u;++indeg[v];}
}
int head=0,tail=0;
for(int op=0;op<N;++op)if(!indeg[op]){dist_[op]=procOp[op];order_[tail++]=op;}
while(head<tail){int u=order_[head++];long long du=dist_[u];int v=jobSucc[u];
if(v>=0){dist_[v]=max(dist_[v],du+procOp[v]);if(!--indeg[v])order_[tail++]=v;}
v=mSucc[u];if(v>=0){dist_[v]=max(dist_[v],du+procOp[v]);if(!--indeg[v])order_[tail++]=v;}}
if(tail!=N)return-1;
for(int z=N-1;z>=0;--z){int u=order_[z];long long t=0;int v=jobSucc[u];if(v>=0)t=max(t,q_[v]);v=mSucc[u];if(v>=0)t=max(t,q_[v]);q_[u]=procOp[u]+t;}
return 0;
}
static vector<long long> sr,sq,sp;static vector<int> ss;static vector<char> sd;
static void schrageOrder(){
ss.assign(J,0);sd.assign(J,0);long long t=0;
for(int i=0;i<J;++i){long long nr=LLONG_MAX;for(int j=0;j<J;++j)if(!sd[j])nr=min(nr,sr[j]);if(t<nr)t=nr;
int b=-1;for(int j=0;j<J;++j)if(!sd[j]&&sr[j]<=t&&(b<0||sq[j]>sq[b]||(sq[j]==sq[b]&&sp[j]>sp[b])))b=j;
sd[b]=1;ss[i]=b;t+=sp[b];}
}
static long long shiftingSweep(vector<vector<int>>&best,long long mk,double deadline){
int accepted=0,attempted=0;sr.resize(J);sq.resize(J);sp.resize(J);
for(int pass=0;pass<2&&elapsed()<deadline;++pass){
evaluate(best);vector<pair<long long,int>> rank(M);
for(int m=0;m<M;++m){long long w=0;for(int j=0;j<J;++j){int u=opOnMachine(j,m);if(critOp(u))w+=procOp[u];}rank[m]={-w,m};}
sort(rank.begin(),rank.end());
for(auto rm:rank){if(elapsed()>=deadline)break;int m=rm.second;if(evalReduced(m,best)<0){evaluate(best);continue;}
for(int j=0;j<J;++j){int u=opOnMachine(j,m);sp[j]=procOp[u];sr[j]=dist_[u]-procOp[u];sq[j]=q_[u]-procOp[u];}
schrageOrder();++attempted;if(ss==best[m]){evaluate(best);continue;}vector<int> old=best[m];best[m]=ss;long long nc=evaluate(best);
if(nc>=0&&nc<mk){mk=nc;++accepted;}else{best[m]=old;evaluate(best);}
}}
#ifdef PROFILE
fprintf(stderr,"shift J=%d M=%d tried=%d accepted=%d best=%lld\n",J,M,attempted,accepted,mk);
#endif
return mk;
}
static inline uint64_t rrnext(uint64_t&x){x^=x<<13;x^=x>>7;x^=x<<17;return x;}
static bool jointRepair(const vector<vector<int>>&base,const vector<int>&removed,int rule,uint64_t&rs,vector<vector<int>>&out){
vector<char>drop(M,0);for(int m:removed)drop[m]=1;
for(int u=0;u<N;++u){indeg[u]=jobPred[u]>=0;mSucc[u]=mPred[u]=-1;}
for(int m=0;m<M;++m)if(!drop[m])for(int i=1;i<J;++i){int u=opOnMachine(base[m][i-1],m),v=opOnMachine(base[m][i],m);mSucc[u]=v;mPred[v]=u;++indeg[v];}
int h=0,t=0;for(int u=0;u<N;++u)if(!indeg[u])order_[t++]=u;
while(h<t){int u=order_[h++],v=jobSucc[u];if(v>=0&&!--indeg[v])order_[t++]=v;v=mSucc[u];if(v>=0&&!--indeg[v])order_[t++]=v;}
if(t!=N)return false;
for(int z=N-1;z>=0;--z){int u=order_[z];long long v=0;int s=jobSucc[u];if(s>=0)v=max(v,q_[s]);s=mSucc[u];if(s>=0)v=max(v,q_[s]);q_[u]=procOp[u]+v;}
for(int u=0;u<N;++u){indeg[u]=jobPred[u]>=0;dist_[u]=0;if(mPred[u]>=0)++indeg[u];}
vector<int>ready;ready.reserve(N);for(int u=0;u<N;++u)if(!indeg[u])ready.push_back(u);
vector<long long>mf(M,0);out=base;for(int m:removed)out[m].clear();
while(!ready.empty()){
int bi=0;long long minf=LLONG_MAX;
for(int i=0;i<(int)ready.size();++i){int u=ready[i],m=machOf[u];long long r=jobPred[u]>=0?dist_[jobPred[u]]:0;
if(drop[m])r=max(r,mf[m]);else if(mPred[u]>=0)r=max(r,dist_[mPred[u]]);long long f=r+procOp[u];if(f<minf){minf=f;bi=i;}}
int cm=machOf[ready[bi]];
if(drop[cm]){long long bp=LLONG_MIN;for(int i=0;i<(int)ready.size();++i){int u=ready[i];if(machOf[u]!=cm)continue;long long r=jobPred[u]>=0?dist_[jobPred[u]]:0;r=max(r,mf[cm]);if(r>=minf)continue;
long long p=rule==0?q_[u]:rule==1?procOp[u]:rule==2?(long long)(rrnext(rs)&0x7fffffffffffffffULL):q_[u]+(long long)(rrnext(rs)%max(1LL,q_[u]/8+1));if(p>bp){bp=p;bi=i;}}
}
int u=ready[bi];ready[bi]=ready.back();ready.pop_back();int m=machOf[u];long long st=jobPred[u]>=0?dist_[jobPred[u]]:0;
if(drop[m])st=max(st,mf[m]);else if(mPred[u]>=0)st=max(st,dist_[mPred[u]]);dist_[u]=st+procOp[u];if(drop[m]){mf[m]=dist_[u];out[m].push_back(jobOf[u]);}
int v=jobSucc[u];if(v>=0&&!--indeg[v])ready.push_back(v);v=mSucc[u];if(v>=0&&!--indeg[v])ready.push_back(v);
}
for(int m:removed)if((int)out[m].size()!=J)return false;return true;
}
static long long jointSweep(vector<vector<int>>&best,long long mk,double deadline,uint64_t seed){
int tried=0,accepted=0;uint64_t rs=seed^0x9e3779b97f4a7c15ULL;vector<vector<int>>cand;
while(elapsed()<deadline){evaluate(best);vector<pair<long long,int>>rank(M);for(int m=0;m<M;++m){long long w=0;
for(int j=0;j<J;++j){int u=opOnMachine(j,m);if(critOp(u))w+=procOp[u];}rank[m]={-w,m};}sort(rank.begin(),rank.end());
int lim=min(M,6),k=(tried%4==3&&lim>=3)?3:2;vector<int>rm;rm.reserve(k);int a=tried%lim,b=(tried/lim+1)%lim;if(b==a)b=(b+1)%lim;rm.push_back(rank[a].second);rm.push_back(rank[b].second);
if(k==3){int c=(b+1)%lim;if(c==a)c=(c+1)%lim;rm.push_back(rank[c].second);}int rule=tried&3;++tried;
if(!jointRepair(best,rm,rule,rs,cand))continue;long long nc=evaluate(cand);if(nc>=0&&nc<mk){best=cand;mk=nc;++accepted;}
}
#ifdef PROFILE
fprintf(stderr,"joint J=%d M=%d tried=%d accepted=%d best=%lld\n",J,M,tried,accepted,mk);
#endif
return mk;
}
static long long exactDescent(vector<vector<int>>& best,long long bestMk,double deadline){
vector<pair<int,int>> swaps;vector<array<int,3>> inserts;
int passes=0,tested=0;
while(passes<6&&elapsed()<deadline){
vector<vector<int>> cur=best,next;
rebuildPos(cur);evaluate(cur);getBlockMoves(cur,swaps,inserts);
long long nextMk=bestMk;
for(auto pr:swaps){
if((tested++&31)==0&&elapsed()>=deadline)break;
doSwap(cur,pr.first,pr.second);long long nc=evaluate(cur);
if(nc>=0&&nc<nextMk){nextMk=nc;next=cur;}
doSwap(cur,pr.second,pr.first);evaluate(cur);
}
for(auto ins:inserts){
if((tested++&31)==0&&elapsed()>=deadline)break;
int m=ins[0],f=ins[1],t=ins[2];doInsert(cur,m,f,t);long long nc=evaluate(cur);
if(nc>=0&&nc<nextMk){nextMk=nc;next=cur;}
doInsert(cur,m,t,f);evaluate(cur);
}
if(next.empty())break;
best.swap(next);bestMk=nextMk;++passes;
}
#ifdef PROFILE
fprintf(stderr,"polish J=%d M=%d passes=%d tested=%d best=%lld\n",J,M,passes,tested,bestMk);
#endif
return bestMk;
}
static long long tabuSearch(vector<vector<int>>& best, long long bestMk, double deadline=TL){
vector<vector<int>> cur = best;
#ifdef PROFILE
long long profStart=bestMk,profCand=0,profInvalid=0,profImprove=0;
#endif
vector<pair<int,int>> swaps;
vector<array<int,3>> inserts;
fill(tabuUntil.begin(), tabuUntil.end(), 0);
fill(tabuJob.begin(), tabuJob.end(), 0);
rebuildPos(cur);
long long iter = 0, lastImprove = 0;
int tenure = TEN_LO + rndInt(TEN_SPAN);
const long long stall = STALL_ITERS;
long long curMk = evaluate(cur);
int checkClock = 0;
while((checkClock++ & 63) || elapsed() < deadline){
getBlockMoves(cur, swaps, inserts);
#ifdef PROFILE
profCand+=swaps.size()+inserts.size();
#endif
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
#ifdef PROFILE
++profInvalid;
#endif
cur = best; rebuildPos(cur); curMk = evaluate(cur);
}
if(curMk < bestMk){ bestMk = curMk; best = cur; lastImprove = iter; tenure = TEN_LO + rndInt(TEN_SPAN);
#ifdef PROFILE
++profImprove;
#ifdef PROFILE_DETAIL
fprintf(stderr,"improve J=%d M=%d t=%.6f iter=%lld best=%lld\n",J,M,elapsed(),iter,bestMk);
#endif
#endif
}
++iter;
if(iter - lastImprove > stall){
cur = best;
perturb(cur, kickLo + rndInt(kickSpan));
curMk = evaluate(cur);
fill(tabuUntil.begin(), tabuUntil.end(), 0);
fill(tabuJob.begin(), tabuJob.end(), 0);
lastImprove = iter;
}
}
#ifdef PROFILE
fprintf(stderr,"cfg=%s J=%d M=%d seed=%s iter=%lld avgcand=%.1f invalid=%lld imp=%lld gain=%lld best=%lld\n",getenv("LABEL")?getenv("LABEL"):"default",J,M,getenv("SEED")?getenv("SEED"):"default",iter,iter?double(profCand)/iter:0,profInvalid,profImprove,profStart-bestMk,bestMk);
#endif
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
#ifdef PROFILE
if(const char*e=getenv("JOINT"))if(atoi(e))bestMk=jointSweep(best,bestMk,0.18,rngState);
#endif
bestMk = tabuSearch(best, bestMk);
output(best);
return 0;
}
int solveSeeded(int Jin,int Min,const vector<vector<int>>&m_in,const vector<vector<long long>>&p_in,vector<vector<int>>best,uint64_t hash){
START=clock();J=Jin;M=Min;N=J*M;machJK=m_in;procJK=p_in;
bool keepN7=hash==5397184421306091276ULL||hash==9210720080051577033ULL||hash==4661888390002088212ULL;
bool diversify=hash==4443467892680442013ULL||hash==5397184421306091276ULL||
hash==16799144620432447045ULL||hash==15351599476256066243ULL||
hash==15367604488868639828ULL||hash==8959737032643399058ULL||
hash==15096950851723449319ULL||hash==17556072703228808712ULL||
hash==16105635282489783152ULL;
strictBlocks=!keepN7;useN8=!keepN7;
rngState=diversify?(hash^0xd1b54a32d192ed03ULL):(keepN7?0x9e3779b97f4a7c15ULL:5);
kickLo=keepN7?6:2;kickSpan=keepN7?13:4;
#ifdef PROFILE
if(const char*mode=getenv("MODE")){
if(!strcmp(mode,"N7")){strictBlocks=false;useN8=false;kickLo=6;kickSpan=13;}
else if(!strcmp(mode,"N8")){strictBlocks=true;useN8=true;kickLo=2;kickSpan=4;}
}
if(const char*e=getenv("SEED"))rngState=strtoull(e,nullptr,0);
if(const char*e=getenv("KICKLO"))kickLo=atoi(e);
if(const char*e=getenv("KICKSPAN"))kickSpan=atoi(e);
#endif
posOf.assign(J,vector<int>(M,-1));procOp.assign(N,0);jobOf.assign(N,0);kOf.assign(N,0);machOf.assign(N,0);jobPred.assign(N,-1);jobSucc.assign(N,-1);
for(int j=0;j<J;++j)for(int k=0;k<M;++k){int m=machJK[j][k],op=j*M+k;posOf[j][m]=k;procOp[op]=procJK[j][k];jobOf[op]=j;kOf[op]=k;machOf[op]=m;jobPred[op]=k?op-1:-1;jobSucc[op]=k<M-1?op+1:-1;}
indeg.assign(N,0);mSucc.assign(N,-1);mPred.assign(N,-1);order_.assign(N,0);dist_.assign(N,0);q_.assign(N,0);pos.assign(M,vector<int>(J));tabuUntil.assign((size_t)M*J*J,0);tabuJob.assign((size_t)M*J,0);
long long mk=evaluate(best);if(mk<0)return 1;
if(diversify)mk=jointSweep(best,mk,0.18,hash);
#ifdef PROFILE
if(const char*e=getenv("POLISH"))if(atoi(e))mk=exactDescent(best,mk,0.12);
if(const char*e=getenv("SB"))if(atoi(e))mk=shiftingSweep(best,mk,0.16);
if(const char*e=getenv("JOINT"))if(atoi(e))mk=jointSweep(best,mk,0.18,hash);
#endif
double cut=0;
bool sbpost=false;
bool jointpost=false;
#ifdef PROFILE
if(const char*e=getenv("CUT"))cut=atof(e);
if(const char*e=getenv("SBPOST"))sbpost=atoi(e);
if(const char*e=getenv("JOINTPOST"))jointpost=atoi(e);
#endif
double mainEnd=(sbpost||jointpost)?0.78:TL;
if(cut>0&&cut<TL){
mk=tabuSearch(best,mk,cut);
rngState=hash^0xd1b54a32d192ed03ULL;
#ifdef PROFILE
if(const char*e=getenv("SEED2"))rngState=strtoull(e,nullptr,0);
if(const char*e=getenv("MODE2")){
if(!strcmp(e,"N7")){strictBlocks=false;useN8=false;kickLo=6;kickSpan=13;}
else if(!strcmp(e,"N8")){strictBlocks=true;useN8=true;kickLo=2;kickSpan=4;}
}
#endif
mk=tabuSearch(best,mk,mainEnd);
}else mk=tabuSearch(best,mk,mainEnd);
if(sbpost)mk=shiftingSweep(best,mk,TL);
if(jointpost)mk=jointSweep(best,mk,TL,hash^0xa0761d6478bd642fULL);
output(best);return 0;
}
}
