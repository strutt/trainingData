#include "BaseList.h"

static std::vector<BaseList::base> baseList;
static bool initializedBaseList = false;

void makeBaseList(){
  baseList.clear();
  baseList.reserve(150); // approx
  baseList.push_back(BaseList::base("Belgrano II", "COMNAP",-77.8744444,-34.6269444,250));
  baseList.push_back(BaseList::base("Brown", "COMNAP",-64.89536982,-62.87045109,10));
  baseList.push_back(BaseList::base("Cámara", "COMNAP",-62.59385833,-59.91934444,22));
  baseList.push_back(BaseList::base("Decepcíon", "COMNAP",-62.97676018,-60.70069463,7));
  baseList.push_back(BaseList::base("Esperanza", "COMNAP",-63.39695938,-56.99805434,25));
  baseList.push_back(BaseList::base("Carlini (formally known as Jubany)", "COMNAP",-62.2379,-58.66684167,10));
  baseList.push_back(BaseList::base("Marambio", "COMNAP",-64.24176895,-56.62322459,200));
  baseList.push_back(BaseList::base("Matienzo", "COMNAP",-64.97586484,-60.07094839,32));
  baseList.push_back(BaseList::base("Melchior", "COMNAP",-64.3257045,-62.97632863));
  baseList.push_back(BaseList::base("Orcadas", "COMNAP",-60.7375919,-44.73738674,4));
  baseList.push_back(BaseList::base("Petrel", "COMNAP",-63.47830316,-56.23099341,18));
  baseList.push_back(BaseList::base("Primavera", "COMNAP",-64.1558536,-60.95425641,50));
  baseList.push_back(BaseList::base("San Martín", "COMNAP",-68.13030011,-67.10293129,5));
  baseList.push_back(BaseList::base("Edgeworth David", "COMNAP",-66.24993451,100.6042134,15));
  baseList.push_back(BaseList::base("Wilkins Aerodrome", "COMNAP",-66.68966697,111.4845668,740));
  baseList.push_back(BaseList::base("Casey", "COMNAP",-66.28234266,110.5267924,30));
  baseList.push_back(BaseList::base("Davis", "COMNAP",-68.57593768,77.96951603,15));
  baseList.push_back(BaseList::base("Mawson", "COMNAP",-67.60264444,62.87302778,5));
  baseList.push_back(BaseList::base("Beaver Lake", "COMNAP",-70.803,68.17983333));
  baseList.push_back(BaseList::base("Law - Racovita - Negoita", "COMNAP",-69.3882965,76.38069645,65));
  baseList.push_back(BaseList::base("Princess Elisabeth", "COMNAP",-71.9498578,23.34689116));
  baseList.push_back(BaseList::base("Comandante Ferraz", "COMNAP",-62.08462076,-58.39256957,8));
  baseList.push_back(BaseList::base("Ohridski", "COMNAP",-62.64072307,-60.36523384,10));
  baseList.push_back(BaseList::base("Lieutenant Arturo Parodi", "COMNAP",-80.31193304,-81.36657271,880));
  baseList.push_back(BaseList::base("Lieutenant Rodolfo Marsh M. Aerodrome", "COMNAP", -62.19373488,-58.9799737,45));
  baseList.push_back(BaseList::base("Arturo Prat", "COMNAP",-62.47933555,-59.66351347,10));
  baseList.push_back(BaseList::base("Lieutenant Luis Carvajal Villarroel", "COMNAP",-67.76132249,-68.91481561));
  baseList.push_back(BaseList::base("Julio Escudero", "COMNAP",-62.20137307,-58.96267648,10));
  baseList.push_back(BaseList::base("Eduardo Frei Montalva", "COMNAP",-62.20023174,-58.96262712,10));
  baseList.push_back(BaseList::base("Bernardo O'Higgins Riquelme", "COMNAP",-63.32095118,-57.89978113,12));
  baseList.push_back(BaseList::base("Ripamonti", "COMNAP",-62.21015169,-58.93474027,50));
  baseList.push_back(BaseList::base("Risopatrón", "COMNAP",-62.37853344,-59.70072429,40));
  baseList.push_back(BaseList::base("President Gabriel Gonzalez Videla", "COMNAP",-64.82386172,-62.85749916));
  baseList.push_back(BaseList::base("Guillermo Mann", "COMNAP",-62.46668333,-60.76668333));
  baseList.push_back(BaseList::base("Sub Base Yelcho", "COMNAP",-64.87585799,-63.58380638,10));
  baseList.push_back(BaseList::base("Great Wall", "COMNAP",-62.21638889,-58.9644444,10));
  baseList.push_back(BaseList::base("Kunlun", "COMNAP",-80.41694444,77.1161111,4087));
  baseList.push_back(BaseList::base("Taishan", "COMNAP",-73.86388333,76.97471666,2620));
  baseList.push_back(BaseList::base("Zhongshan", "COMNAP",-69.37333333,76.37777778,10));
  baseList.push_back(BaseList::base("Johann Gregor Mendel", "COMNAP",-63.80062531,-57.8825918,10));
  baseList.push_back(BaseList::base("Refugio Ecuador", "COMNAP",-62.12099018,-58.3951612,10));
  baseList.push_back(BaseList::base("Maldonado", "COMNAP",-62.44932928,-59.74096859,10));
  baseList.push_back(BaseList::base("Aboa", "COMNAP",-73.04228194,-13.40735764,400));
  baseList.push_back(BaseList::base("Cap Prud'homme", "COMNAP",-66.68760751,139.9071699,10));
  baseList.push_back(BaseList::base("Dumont d'Urville", "COMNAP",-66.66283333,140.0013333,42));
  baseList.push_back(BaseList::base("Concordia", "COMNAP",-75.0999745,123.3326053,3220));
  baseList.push_back(BaseList::base("Dallman Lab at Base Carlini", "COMNAP",-62.23760854,-58.66671638));
  baseList.push_back(BaseList::base("Antarctic Receiving Station (GARS)", "COMNAP",-63.32109182,-57.90094008));
  baseList.push_back(BaseList::base("Gondwana", "COMNAP",-74.63550051,164.2212145));
  baseList.push_back(BaseList::base("Kohnen", "COMNAP",-75.00191485,0.066733305,2900));
  baseList.push_back(BaseList::base("Neumayer III", "COMNAP",-70.67725032,-8.271599716,40));
  baseList.push_back(BaseList::base("Bharati", "COMNAP",-69.4,76.183));
  baseList.push_back(BaseList::base("Dakshin Gangotri", "COMNAP",-70.08333333,12.0));
  baseList.push_back(BaseList::base("Maitri", "COMNAP",-70.76683367,11.73078318,130));
  baseList.push_back(BaseList::base("Browning Pass", "COMNAP",-74.62290575,163.9151708,170));
  baseList.push_back(BaseList::base("Enigma Lake", "COMNAP",-74.719164,164.027711,170));
  baseList.push_back(BaseList::base("Mid Point", "COMNAP",-75.54170266,145.8203993,2520));
  baseList.push_back(BaseList::base("Sitry", "COMNAP",-71.65173817,148.654737,1600));
  baseList.push_back(BaseList::base("Mario Zucchelli", "COMNAP",-74.69480668,164.1132678,15));
  baseList.push_back(BaseList::base("S17", "COMNAP",-69.02807608,40.09244947,620));
  baseList.push_back(BaseList::base("Asuka", "COMNAP",-71.52627263,24.11251621));
  baseList.push_back(BaseList::base("Dome Fuji", "COMNAP",-77.31707524,39.69877019,3810));
  baseList.push_back(BaseList::base("Mizuho", "COMNAP",-70.6990456,44.27781168));
  baseList.push_back(BaseList::base("Syowa", "COMNAP",-69.00412231,39.58183616,29));
  baseList.push_back(BaseList::base("King Sejong", "COMNAP",-62.2232395,-58.7865044,10));
  baseList.push_back(BaseList::base("Jang Bogo", "COMNAP",-74.61666667,164.2));
  baseList.push_back(BaseList::base("Dirck Gerritsz Laboratory ", "COMNAP",-67.56863642,-68.124377,16));
  baseList.push_back(BaseList::base("Scott Base", "COMNAP",-77.849437,166.767281,10));
  baseList.push_back(BaseList::base("Tor", "COMNAP",-71.88951483,5.159892622,1625));
  baseList.push_back(BaseList::base("Troll", "COMNAP",-72.01194552,2.533070031,1300));
  baseList.push_back(BaseList::base("Machu Picchu", "COMNAP",-62.09159252,-58.47055896,10));
  baseList.push_back(BaseList::base("Arctowski", "COMNAP",-62.15977292,-58.47332175,2));
  baseList.push_back(BaseList::base("Molodezhnaya Airfield", "COMNAP",-67.66974075,45.82860961,225));
  baseList.push_back(BaseList::base("Novolazarevskaya Airfield", "COMNAP",-70.8217966,11.63987117,550));
  baseList.push_back(BaseList::base("Bellingshausen", "COMNAP",-62.19818377,-58.96060891,16));
  baseList.push_back(BaseList::base("Druzhnaya-4", "COMNAP",-69.74782719,73.70921376,20));
  baseList.push_back(BaseList::base("Leningradskaya", "COMNAP",-69.50149558,159.3911476));
  baseList.push_back(BaseList::base("Mirny", "COMNAP",-66.55823444,93.0003795,40));
  baseList.push_back(BaseList::base("Molodezhnaya", "COMNAP",-67.66544285,45.84203764,42));
  baseList.push_back(BaseList::base("Novolazarevskaya", "COMNAP",-70.77692572,11.82366729,102));
  baseList.push_back(BaseList::base("Progress", "COMNAP",-69.37812036,76.38779398,15));
  baseList.push_back(BaseList::base("Russkaya", "COMNAP",-74.76572579,-136.8000722));
  baseList.push_back(BaseList::base("Soyuz", "COMNAP",-70.5765771,68.79493694,336));
  baseList.push_back(BaseList::base("Vostok", "COMNAP",-78.46415952,106.8379632,3500));
  baseList.push_back(BaseList::base("SANAE IV", "COMNAP",-71.67287313,-2.840324047,850));
  baseList.push_back(BaseList::base("Gabriel de Castilla", "COMNAP",-62.97719656,-60.67574545,15));
  baseList.push_back(BaseList::base("Juan Carlos I", "COMNAP",-62.66340863,-60.38814603,12));
  baseList.push_back(BaseList::base("Svea", "COMNAP",-74.57606102,-11.22487242));
  baseList.push_back(BaseList::base("Wasa", "COMNAP",-73.04279886,-13.41290971,400));
  baseList.push_back(BaseList::base("Vernadsky", "COMNAP",-65.24574295,-64.25748123,7));
  baseList.push_back(BaseList::base("Fossil Bluff", "COMNAP",-71.32341152,-68.28906406,92));
  baseList.push_back(BaseList::base("Rothera Skiway", "COMNAP",-67.56758573,-68.12717991,250));
  baseList.push_back(BaseList::base("Sky Blu", "COMNAP",-74.856351,-71.58519987,1370-1500));
  baseList.push_back(BaseList::base("Halley", "COMNAP",-75.57982394,-26.72862432,37));
  baseList.push_back(BaseList::base("Rothera", "COMNAP",-67.56863642,-68.124377,16));
  baseList.push_back(BaseList::base("Signy", "COMNAP",-60.70828739,-45.59538916,5));
  baseList.push_back(BaseList::base("Artigas", "COMNAP",-62.18455071,-58.90244173,17));
  baseList.push_back(BaseList::base("Ruperto Elichiribehety", "COMNAP",-63.40237211,-56.9909068));
  baseList.push_back(BaseList::base("Marble Point Heliport", "COMNAP",-77.41366667,163.6791667));
  baseList.push_back(BaseList::base("Odell Glacier Camp", "COMNAP",-76.6601162,159.953156,1600));
  baseList.push_back(BaseList::base("Siple Dome", "COMNAP",-81.65430285,-149.0051336));
  baseList.push_back(BaseList::base("Amundsen-Scott South Pole Station", "COMNAP",-89.9975,139.2728,2830));
  baseList.push_back(BaseList::base("McMurdo Station", "COMNAP",-77.8482091,166.668422,10));
  baseList.push_back(BaseList::base("Palmer Station", "COMNAP",-64.77426073,-64.05333352,10));





  baseList.push_back(BaseList::base("WAIS divide", "http://waisdivide.unh.edu/about/site.shtml#selection",
				    -79.467, -112.085, 1766));


  initializedBaseList = true;
}

const BaseList::base& BaseList::getBase(UInt_t index){
  if(!initializedBaseList){
    makeBaseList();
    initializedBaseList = true;
  }

  index = index < baseList.size() ? index : 0;
  return baseList.at(index);
}


size_t BaseList::getNumBases(){
  if(!initializedBaseList){
    makeBaseList();
    initializedBaseList = true;
  }
  return baseList.size();
}
