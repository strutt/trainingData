#include "TString.h"
#include <vector>

bool initializedBaseList = false;

namespace COMNAP2014{
  // a list of all the bases


  class base{
  public:
    base(const TString& theName, double theLat, double theLon, double theAlt=0){
      name = theName;
      lat = theLat;
      lon = theLon;
      alt = theAlt;
    }
    TString name;
    double lat;
    double lon;
    double alt;
  };

  std::vector<base> baseList;
  // std::map<TString, base> baseList;


  void makeBaseList(){
    baseList.clear();
    baseList.reserve(150); // approx
    baseList.push_back(base("Belgrano II",-77.8744444,-34.6269444,250));
    baseList.push_back(base("Brown",-64.89536982,-62.87045109,10));
    baseList.push_back(base("Cámara",-62.59385833,-59.91934444,22));
    baseList.push_back(base("Decepcíon",-62.97676018,-60.70069463,7));
    baseList.push_back(base("Esperanza",-63.39695938,-56.99805434,25));
    baseList.push_back(base("Carlini (formally known as Jubany)",-62.2379,-58.66684167,10));
    baseList.push_back(base("Marambio",-64.24176895,-56.62322459,200));
    baseList.push_back(base("Matienzo",-64.97586484,-60.07094839,32));
    baseList.push_back(base("Melchior",-64.3257045,-62.97632863));
    baseList.push_back(base("Orcadas",-60.7375919,-44.73738674,4));
    baseList.push_back(base("Petrel",-63.47830316,-56.23099341,18));
    baseList.push_back(base("Primavera",-64.1558536,-60.95425641,50));
    baseList.push_back(base("San Martín",-68.13030011,-67.10293129,5));
    baseList.push_back(base("Edgeworth David",-66.24993451,100.6042134,15));
    baseList.push_back(base("Wilkins Aerodrome",-66.68966697,111.4845668,740));
    baseList.push_back(base("Casey",-66.28234266,110.5267924,30));
    baseList.push_back(base("Davis",-68.57593768,77.96951603,15));
    baseList.push_back(base("Mawson",-67.60264444,62.87302778,5));
    baseList.push_back(base("Beaver Lake",-70.803,68.17983333));
    baseList.push_back(base("Law - Racovita - Negoita",-69.3882965,76.38069645,65));
    baseList.push_back(base("Princess Elisabeth",-71.9498578,23.34689116));
    baseList.push_back(base("Comandante Ferraz",-62.08462076,-58.39256957,8));
    baseList.push_back(base("Ohridski",-62.64072307,-60.36523384,10));
    baseList.push_back(base("Lieutenant Arturo Parodi",-80.31193304,-81.36657271,880));
    baseList.push_back(base("Lieutenant Rodolfo Marsh M. Aerodrome", -62.19373488,-58.9799737,45));
    baseList.push_back(base("Arturo Prat",-62.47933555,-59.66351347,10));
    baseList.push_back(base("Lieutenant Luis Carvajal Villarroel",-67.76132249,-68.91481561));
    baseList.push_back(base("Julio Escudero",-62.20137307,-58.96267648,10));
    baseList.push_back(base("Eduardo Frei Montalva",-62.20023174,-58.96262712,10));
    baseList.push_back(base("Bernardo O'Higgins Riquelme",-63.32095118,-57.89978113,12));
    baseList.push_back(base("Ripamonti",-62.21015169,-58.93474027,50));
    baseList.push_back(base("Risopatrón",-62.37853344,-59.70072429,40));
    baseList.push_back(base("President Gabriel Gonzalez Videla",-64.82386172,-62.85749916));
    baseList.push_back(base("Guillermo Mann",-62.46668333,-60.76668333));
    baseList.push_back(base("Sub Base Yelcho",-64.87585799,-63.58380638,10));
    baseList.push_back(base("Great Wall",-62.21638889,-58.9644444,10));
    baseList.push_back(base("Kunlun",-80.41694444,77.1161111,4087));
    baseList.push_back(base("Taishan",-73.86388333,76.97471666,2620));
    baseList.push_back(base("Zhongshan",-69.37333333,76.37777778,10));
    baseList.push_back(base("Johann Gregor Mendel",-63.80062531,-57.8825918,10));
    baseList.push_back(base("Refugio Ecuador",-62.12099018,-58.3951612,10));
    baseList.push_back(base("Maldonado",-62.44932928,-59.74096859,10));
    baseList.push_back(base("Aboa",-73.04228194,-13.40735764,400));
    baseList.push_back(base("Cap Prud'homme",-66.68760751,139.9071699,10));
    baseList.push_back(base("Dumont d'Urville",-66.66283333,140.0013333,42));
    baseList.push_back(base("Concordia",-75.0999745,123.3326053,3220));
    baseList.push_back(base("Dallman Lab at Base Carlini",-62.23760854,-58.66671638));
    baseList.push_back(base("Antarctic Receiving Station (GARS)",-63.32109182,-57.90094008));
    baseList.push_back(base("Gondwana",-74.63550051,164.2212145));
    baseList.push_back(base("Kohnen",-75.00191485,0.066733305,2900));
    baseList.push_back(base("Neumayer III",-70.67725032,-8.271599716,40));
    baseList.push_back(base("Bharati",-69.4,76.183));
    baseList.push_back(base("Dakshin Gangotri",-70.08333333,12.0));
    baseList.push_back(base("Maitri",-70.76683367,11.73078318,130));
    baseList.push_back(base("Browning Pass",-74.62290575,163.9151708,170));
    baseList.push_back(base("Enigma Lake",-74.719164,164.027711,170));
    baseList.push_back(base("Mid Point",-75.54170266,145.8203993,2520));
    baseList.push_back(base("Sitry",-71.65173817,148.654737,1600));
    baseList.push_back(base("Mario Zucchelli",-74.69480668,164.1132678,15));
    baseList.push_back(base("S17",-69.02807608,40.09244947,620));
    baseList.push_back(base("Asuka",-71.52627263,24.11251621));
    baseList.push_back(base("Dome Fuji",-77.31707524,39.69877019,3810));
    baseList.push_back(base("Mizuho",-70.6990456,44.27781168));
    baseList.push_back(base("Syowa",-69.00412231,39.58183616,29));
    baseList.push_back(base("King Sejong",-62.2232395,-58.7865044,10));
    baseList.push_back(base("Jang Bogo",-74.61666667,164.2));
    baseList.push_back(base("Dirck Gerritsz Laboratory ",-67.56863642,-68.124377,16));
    baseList.push_back(base("Scott Base",-77.849437,166.767281,10));
    baseList.push_back(base("Tor",-71.88951483,5.159892622,1625));
    baseList.push_back(base("Troll",-72.01194552,2.533070031,1300));
    baseList.push_back(base("Machu Picchu",-62.09159252,-58.47055896,10));
    baseList.push_back(base("Arctowski",-62.15977292,-58.47332175,2));
    baseList.push_back(base("Molodezhnaya Airfield",-67.66974075,45.82860961,225));
    baseList.push_back(base("Novolazarevskaya Airfield",-70.8217966,11.63987117,550));
    baseList.push_back(base("Bellingshausen",-62.19818377,-58.96060891,16));
    baseList.push_back(base("Druzhnaya-4",-69.74782719,73.70921376,20));
    baseList.push_back(base("Leningradskaya",-69.50149558,159.3911476));
    baseList.push_back(base("Mirny",-66.55823444,93.0003795,40));
    baseList.push_back(base("Molodezhnaya",-67.66544285,45.84203764,42));
    baseList.push_back(base("Novolazarevskaya",-70.77692572,11.82366729,102));
    baseList.push_back(base("Progress",-69.37812036,76.38779398,15));
    baseList.push_back(base("Russkaya",-74.76572579,-136.8000722));
    baseList.push_back(base("Soyuz",-70.5765771,68.79493694,336));
    baseList.push_back(base("Vostok",-78.46415952,106.8379632,3500));
    baseList.push_back(base("SANAE IV",-71.67287313,-2.840324047,850));
    baseList.push_back(base("Gabriel de Castilla",-62.97719656,-60.67574545,15));
    baseList.push_back(base("Juan Carlos I",-62.66340863,-60.38814603,12));
    baseList.push_back(base("Svea",-74.57606102,-11.22487242));
    baseList.push_back(base("Wasa",-73.04279886,-13.41290971,400));
    baseList.push_back(base("Vernadsky",-65.24574295,-64.25748123,7));
    baseList.push_back(base("Fossil Bluff",-71.32341152,-68.28906406,92));
    baseList.push_back(base("Rothera Skiway",-67.56758573,-68.12717991,250));
    baseList.push_back(base("Sky Blu",-74.856351,-71.58519987,1370-1500));
    baseList.push_back(base("Halley",-75.57982394,-26.72862432,37));
    baseList.push_back(base("Rothera",-67.56863642,-68.124377,16));
    baseList.push_back(base("Signy",-60.70828739,-45.59538916,5));
    baseList.push_back(base("Artigas",-62.18455071,-58.90244173,17));
    baseList.push_back(base("Ruperto Elichiribehety",-63.40237211,-56.9909068));
    baseList.push_back(base("Marble Point Heliport",-77.41366667,163.6791667));
    baseList.push_back(base("Odell Glacier Camp",-76.6601162,159.953156,1600));
    baseList.push_back(base("Siple Dome",-81.65430285,-149.0051336));
    baseList.push_back(base("Amundsen-Scott South Pole Station",-89.9975,139.2728,2830));
    baseList.push_back(base("McMurdo Station",-77.8482091,166.668422,10));
    baseList.push_back(base("Palmer Station",-64.77426073,-64.05333352,10));
  }

  const base& getBase(UInt_t index){
    if(!initializedBaseList){
      makeBaseList();
      initializedBaseList = true;
    }

    index = index < baseList.size() ? index : 0;
    return baseList.at(index);
  }



};
