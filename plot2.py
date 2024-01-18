from mpmath import *
import matplotlib.pyplot as plt
from main import *

gu = [mpf('-0.0151481405249576644549874401106725955068458665888282827353282533916185891260820969277895378389592188408296358598469546735579548155999468192330605858407284'), mpf('-197.354231219767254900269921928982445811541425768698369441691648970286539556983811668332436514171509693819568446832444119259944555984177277492487087670952'), mpf('196245.248943122648564832295809376441500386333376198764122024752136574275560522835157653949947966329186444989363623386915351256826554048049234616887248608'), mpf('-63318755.5938234611427889050158825786356989958178803564174083495985471621231317713917805152672979523572470474794418664762046770057463597204095344876004619'), mpf('10601978059.4923871723862977476895175316906085159207526943538475252403258021410744016959626199571438184351773758888474695149568846112120084479337397601964'), mpf('-1091293416268.14302061220795074916881392541349838095907098787540967869556887197366898643556821391613115307554089936223213422622904332885273500553175080261'), mpf('76027472231628.9490820424575860016376067252553216297272175184953839106357952948815118324877314551077792614672287169729614221984182367749274139000144640402'), mpf('-3818019501135481.57210830206124900487925467463976243190968327390231318674288262665959620464139864640572011613147143372804296093082795905445638937695402299'), mpf('144519488874627485.053532167648279101313588779200886650938106276007907331266794997311522093695243114845182413293873651984637024438464614723657170600044622'), mpf('-4262665924266187214.67191220319204162102821783955562112089959053676065304041994039275902794050471688027704913365646820157016858021801346865571929305491612'), mpf('100523247486205297866.49842985375396424742268839014279762404234881420613732906242965015528617718977179972041861667588840765223653141665076181348848657488'), mpf('-1934397585093730307236.2707447656433968870505808200546415284925251964596846974755962066015402186049971193679772787414222475513924443275240847807129946775'), mpf('30881942747077341145533.0091985011176389616077187255856407171840761554790929502334163160082528750182507269763313365089420953514657604425055718483579944555'), mpf('-414633463572693263375047.589368256978998092386912549084643338844794126440787025874809203407788628453371573727636935736626563436491549724259925417317350541'), mpf('4735520487328893051169583.48450831690463300262975181192245825536593282117839735345865651806950822609373564205627361295195202021257023411703373212677027117'), mpf('-46449706098568802689041147.6886323708234417554045391422065841154352564663459450772218201666759796430790615047707296801206752652241335660112787369856479667'), mpf('394505461612161799682243468.718262572306762380279584797834817662036741281007117302752667338637820694628315236937027553053602250437251443609589011325676797'), mpf('-2921507040080020861718707793.16740670157612086338880888137804762957128591664977485133438926739568472656886824662707790279257735199837730977948293570518043'), mpf('18977748291438453231507818399.2161961228513304072887148852431257876887877036236602065068789047281218097546360432384178779646445001057780294339421217228039'), mpf('-108693212318236836670442072907.742031466534618216651728575701716194382353005772959055506877423806039467028951763620876943636551240584520795559590226307074'), mpf('551319584408211272337033793264.979309525998479564185433297647979537443528955774668473036679415348639325979575854309713321487604681332439279659382176085797'), mpf('-2485981714618060370929370084880.74100520980206898270070975880414327301049185028130515768031695864104723447627517439637612010801343783345364113915730756424'), mpf('9997573965248854064815110432226.96123229612431494511082707043806723192738142465932351592431119617194762493582193491617877746584730863034310424133047727608'), mpf('-35957149885683703558223169829686.2002488155821135613636591893790520610713779762485926061527805538287686364984637138058971846552251137091725078941560979486'), mpf('115921464623886784869031905481446.980461764567146903664325674894786597137574905420176414263125015052920951140774803565483466194332103988583595028813140885'), mpf('-335614992039545343164764388303240.609073636219825823240303161310886102044006131062438802479432012550806481960054361650091703485012836183080888669139572187'), mpf('873889007519829329754953992258196.219002184622642846356152513808066527436210802100283072687642823792296179483070501046524474888810993058993031490667999123'), mpf('-2048717050209182900657277331634704.6281916043807160805236358276799860272706230525919236611869680073243175073267010375407187933351479102517062561144648161'), mpf('4327435353730488034244036304597836.10095261854633920176333596331149362802957406314014058178436058337984563592893514007798717186738422283483762678389948243'), mpf('-8238606875964043204298527154202485.14629500938167602674355570044767679494625079631653000839062900537741253078393254319352233122151986447040676703156996933'), mpf('14136450040829230578977240955169864.2225387155527440540267289661454192903174084819939230167116949376096646730679030806586422986127005138070824500986361667'), mpf('-21852964645299996917768318137681340.5893508084209896145670321519815500824400427457586177588839448380493760926356378320823543434247186505423076970565713539'), mpf('30409215618656706318717522006015683.39306423160187073289326698320059002651008817226133336119144039287466771894908414483069903820897268609950947022426563'), mpf('-38042879370704897872677869817715028.88100827034112822265733598693343298831825533862254876527393619293930963717193375775618936111577614348148653892169685'), mpf('42712040775669407244847861834875622.403598725463663395506053826464129054321354318394598286806710061624991467737791285336200041969462934094451072158857597'), mpf('-42936976538469986626831642158540520.2415739631841187497717839908843510204797035323232035340542741157172438552307163655938047112403033341292288122064913945'), mpf('38533267292511452110900089265430152.2838842079723849814703303255229473804794580216290493820424110644116234792259296379249534874200849489132117525118978823'), mpf('-30757951085997358161153239713151596.4679918348681119447497280623276573982158931271062457250664055936044980139271634341446117355534461804983765636645050092'), mpf('21736978086760191181003826970601447.1426469477864703294607304926407792753864629014313169017154414773318170655427184249758193805028006788002447621247421727'), mpf('-13523307896130848995964113467078079.9589129985494674932654697644436686317723552559996737976055816898611071839440661691946334731229049156493412946954311548'), mpf('7353992685921960855352867674359790.56497988648957182450007186999822962754979891915259503046578630132442291011215152766697861372746660118806518887258453432'), mpf('-3464546266905257142188755195801119.85792165826647750846934351940752549771373307412522105331720687231849200505462792181948895566834179249989536238019121273'), mpf('1398051973855315551977286727091645.76428290619735456202972728176055183382599528058353639479748217531167674935869934669357498984638815821501849832522314699'), mpf('-476164525752550465327756612427462.818685568572895330764120168532329877372621930978001320991956452066726230669479316289435527533282531547005090195115871321'), mpf('134218781020827434645423493692580.418168970329818246509803873283719357316027034813964076155158862659980755828420885666309443695515216525166103255936821437'), mpf('-30468981785514658275112619830292.8416783120713760095565779696618413791274600104371201866479661478747874213801576081499364683539186295311047509361131949474'), mpf('5352339425538240613795856564797.96513569524133499996154153332334163944275756177715903233634693867557804998506834721320544166935019812315731946991882143914'), mpf('-682661599505646523954878271371.296233485523668575225066437812917903972239736169640833788681413971662875648242394439204977778785361174241613430197899607485'), mpf('56239080539447770118548186004.4693176465112340172928412767863030337318727859685040525384201415689214798882997926504278135938386874362296829785535140036089'), mpf('-2246130660085715581433144968.45251751946297470359818192717321008897706703228304162754862438414148126103698914102423138297560953154643740150268200106441156')]
fu = [mpf('0.999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999389'), mpf('2.61065232906792236661256570307813505976418305184559867571047740157011013097407484792469046710547474454538410837714438038686487398108183172187829788300815e-152'), mpf('-0.166666833477526735442401107435674172220110793387752493778573737037476377227156229328306083673104952596605171462049232302443045378374624175340686060403454'), mpf('0.00000889441236044042188267639332120395449263900212080695654538162389374171395926562092002535006749365199956551735634396548338579940545472764057514177211434529'), mpf('0.0249023997798559703726219366239817959185265679463177691917590522289380299597495753812894306045412760605156512184197231030819219195098962440939405316917044'), mpf('0.000485459934887731876044210464866324508259131040922493201103202699411337014416415174976655146193430476481264162034480184508794086932961123863699641569506732'), mpf('-0.00516812449372417238925614678590691393788412051775365263172347362310617158595863986226526749114720482889396935636947803405936553964376958578099510449372207'), mpf('0.00262110850434118516206180009819249701567671700377308095790973638762197181052257879132064567640382822772124226018114567678425663928235223351710957445335137'), mpf('-0.00284184001539703150082236441347825460087652462044968784005833204338465212257584675607968249138812649629430096302838639019797814283032779764222330108597105'), mpf('0.00319499171997960913311491945766396552683432378279041260190814808171623741292268417377906644249296230557904901569996671302212886986858066586883778234688015'), mpf('-0.00228060951177106618088022919975068597118201465413280986256674868211123420517208259862876350688455105823407344646413953424786890449625251329325377356956165'), mpf('0.00110954012591148180494733230415905935605436086279833873293891218543743801632169373756531906982487766647847516343613648667575917468084103970297151676758548'), mpf('-0.0003958833157830565726340102345636129831601001769028869504588798112059320238590432387352295594331301889805837053355727226840409165575062585308737413177888'), mpf('0.000108217634243560547816480064241692228216105991619558549262267420667894330686179980286644973072355006164626335267064017966015361062946688117112418205128313'), mpf('-0.0000232380021665033159396575150887717718185660210305640870106324866804894587245610183572893354759039253913137314009452711250045143456676730483695402682087965'), mpf('0.00000397329698849560407361290952246403317054183127776071856453681409616582481271613413110269324231015861586951792146387792846131578874545851176730828105423958'), mpf('-0.000000544045643983522431212934559945553576415278060049026411681167160626759702916477812023278070056010464579069221432996142307985053046590513786215179995490004'), mpf('0.0000000596256533182339049461116315378544587275048654256854852809618829498770383607623704119386664128500338273208315041619192762143249103167462949712325001722622'), mpf('-0.00000000519724707037687139963059562416310195821085500975442096075806009008044759323918182644203307846986933714549934518350139037549113882504963951463953793893471'), mpf('0.000000000355632027958458400251764189845353909780586908634500892193814657245746644936682723440548937041524823878488199129781752277563430070922556445287091257238959'), mpf('-0.0000000000186907644408515791426652302004721316042176383788615295715594789497890903451685912265876942000494937246709341198753162826825168272589821555656165464997225'), mpf('0.000000000000728209395921390364917644889123507171307861699102181177784997037773032887392875904809606424076002788832725911927711385035618009792520224715987081155585802'), mpf('-0.0000000000000198096170527680100019030083163969316775774258360350375272866864413149066648068731060621810032328940187487163812272991373297917905109076466423915765910838'), mpf('0.00000000000000033581554172735922852876844794523204649308419234720273468168978719139831487417791091210051224046988544756546353909167820979759729110529684625386532994167'), mpf('-0.00000000000000000267028532997367320311339176702926614913291275929309649065240761613644626205671985808085643061968523524913905799480701004165105941642703934122186138662478')]


plt.clf()
plt.grid()
plt.title("Piecewise Solution")
ax = plt.gca()

ax.set_ylim([-0.5, 1.3])

D = [0,30]
points = 200
dx = (D[1]-D[0])/points
plt.xlabel("x")
plt.ylabel("y")

def f(x):
	if x < 31:
		return polyeval(fu,x)
	else:
		u = exp(-x)
		print(polyeval(gu,u))
		return polyeval(gu,u)



xlist = [D[0]+i*dx for i in range(points+1)]

aylist = [f(x) for x in xlist]
pylist = [polyeval(gu,exp(-x)) for x in xlist]
plt.plot(xlist,pylist,"-r")


plt.plot(xlist,aylist, "-g",label="Piecewise Solution")
plt.plot(xlist,pylist,"--r")


plt.legend()
plt.show()


