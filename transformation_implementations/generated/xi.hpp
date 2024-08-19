/*
 * Generated. DO NOT EDIT.
 */

#ifndef ECEF2GEO_XI_HPP
#define ECEF2GEO_XI_HPP

namespace ecef2geo {
namespace priv {

/** Minimax atan approximations valid on [-1,1].
 *
 * @param N Degree of the polynomial.
 * @param x value in [-1,1].
 * @return An approximation of atan(x).
 */
template<int N> inline double xi(double x) = delete;  // Only allow provided specializations.
template<> inline double xi<1>(double x) { return x * (8.3327886416972108e-01); }
template<> inline double xi<2>(double x) { double x2 = x * x; return x * (9.7239418927023847e-01 + -1.9194804039207914e-01 * x2); }
template<> inline double xi<3>(double x) { double x2 = x * x; double x22 = x2 * x2; return x * (9.9535795470534027e-01 + -2.8869023791774234e-01 * x2 + 7.9339041370619872e-02 * x22); }
template<> inline double xi<4>(double x) { double x2 = x * x; double x22 = x2 * x2; return x * (9.9921381224115369e-01 + -3.2117496590422629e-01 * x2 + (1.4626445659593218e-01 + -3.8986510281795543e-02 * x2) * x22); }
template<> inline double xi<5>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; return x * (9.9986632904710004e-01 + -3.3030478126799267e-01 * x2 + (1.8015928178377459e-01 + -8.5156335459610334e-02 * x2) * x22 + 2.0845107810858079e-02 * x24); }
template<> inline double xi<6>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; return x * (9.9997721903010127e-01 + -3.3262282742336128e-01 * x2 + (1.9354037427200508e-01 + -1.1642647899693130e-01 * x2) * x22 + (5.2647349546537023e-02 + -1.1719135396032095e-02 * x2) * x24); }
template<> inline double xi<7>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; return x * (9.9999611156002193e-01 + -3.3317368072014072e-01 * x2 + (1.9807815660563511e-01 + -1.3233342337191317e-01 * x2) * x22 + (7.9623675398178760e-02 + -3.3604222386222235e-02 * x2 + 6.8117937024772251e-03 * x22) * x24); }
template<> inline double xi<8>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; return x * (9.9999933557895379e-01 + -3.3329860785866898e-01 * x2 + (1.9946565664948636e-01 + -1.3908629610769453e-01 * x2) * x22 + (9.6421974777848821e-02 + -5.5912328807716874e-02 * x2 + (2.1862959296299865e-02 + -4.0545676075818409e-03 * x2) * x22) * x24); }
template<> inline double xi<9>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999988638310751e-01 + -3.3332597028089081e-01 * x2 + (1.9985906766141972e-01 + -1.4161229201159091e-01 * x2) * x22 + (1.0498946125718338e-01 + -7.2348574954356469e-02 * x2 + (3.9781225234399571e-02 + -1.4401358998390296e-02 * x2) * x22) * x24 + 2.4567248628024892e-03 * x28); }
template<> inline double xi<10>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999998056032907e-01 + -3.3333180376928163e-01 * x2 + (1.9996436813152854e-01 + -1.4247222678707639e-01 * x2) * x22 + (1.0878009931861596e-01 + -8.2137617056466841e-02 * x2 + (5.5028101156174181e-02 + -2.8490776169781917e-02 * x2) * x22) * x24 + (9.5673402955380176e-03 + -1.5093031760903022e-03 * x2) * x28); }
template<> inline double xi<11>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999999667245385e-01 + -3.3333302089692613e-01 * x2 + (1.9999129801285570e-01 + -1.4274432171781580e-01 * x2) * x22 + (1.1028651401951308e-01 + -8.7138616815786425e-02 * x2 + (6.5413817590740830e-02 + -4.2088159206687315e-02 * x2) * x22) * x24 + (2.0467886306328867e-02 + -6.3947943012842947e-03 * x2 + 9.3756387415165670e-04 * x22) * x28); }
template<> inline double xi<12>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999999943022189e-01 + -3.3333327040718965e-01 * x2 + (1.9999793539915971e-01 + -1.4282551389892477e-01 * x2) * x22 + (1.1083671099533618e-01 + -8.9411184337301530e-02 * x2 + (7.1430746289771899e-02 + -5.2514556999956749e-02 * x2) * x22) * x24 + (3.2232566609514774e-02 + -1.4721216068835910e-02 * x2 + (4.2936462844677191e-03 + -5.8769992093623553e-04 * x2) * x22) * x28); }
template<> inline double xi<13>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999999990240966e-01 + -3.3333332081059667e-01 * x2 + (1.9999952202232363e-01 + -1.4284860322573315e-01 * x2) * x22 + (1.1102438787732532e-01 + -9.0352255555093739e-02 * x2 + (7.4508238425781743e-02 + -5.9269197198135301e-02 * x2) * x22) * x24 + (4.2261257239378133e-02 + -2.4658752476606615e-02 * x2 + (1.0588543381896344e-02 + -2.8928368099830404e-03 * x2) * x22 + 3.7118062799720641e-04 * x24) * x28); }
template<> inline double xi<14>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999999998328136e-01 + -3.3333333086677668e-01 * x2 + (1.9999989164783083e-01 + -1.4285491033770566e-01 * x2) * x22 + (1.1108488389521914e-01 + -9.0713495063830683e-02 * x2 + (7.5933101852628750e-02 + -6.3109521418777435e-02 * x2) * x22) * x24 + (4.9445357766039380e-02 + -3.3981220148545598e-02 * x2 + (1.8820607710988513e-02 + -7.6115686823698953e-03 * x2) * x22 + (1.9542989485087399e-03 + -2.3593188960412083e-04 * x2) * x24) * x28); }
template<> inline double xi<15>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999999999713534e-01 + -3.3333333285183776e-01 * x2 + (1.9999997588239277e-01 + -1.4285657527852751e-01 * x2) * x22 + (1.1110347572689653e-01 + -9.0843662820425592e-02 * x2 + (7.6540999501698302e-02 + -6.5075415596078731e-02 * x2) * x22) * x24 + (5.3939982359764723e-02 + -4.1298258386912323e-02 * x2 + (2.7248467694973274e-02 + -1.4321564110192034e-02 * x2) * x22 + (5.4663698968691514e-03 + -1.3230853816482112e-03 * x2 + 1.5078676343054553e-04 * x22) * x24) * x28); }
template<> inline double xi<16>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; return x * (9.9999999999950908e-01 + -3.3333333324007675e-01 * x2 + (1.9999999471702140e-01 + -1.4285700204199351e-01 * x2) * x22 + (1.1110896094930681e-01 + -9.0888117301325033e-02 * x2 + (7.6783159568187618e-02 + -6.5998269071639969e-02 * x2) * x22) * x24 + (5.6460742263488186e-02 + -4.6293883569208967e-02 * x2 + (3.4440998399993243e-02 + -2.1766144335727145e-02 * x2) * x22 + (1.0861830271605163e-02 + -3.9212439249656297e-03 * x2 + (8.9729819815698214e-04 + -9.6827484897702223e-05 * x2) * x22) * x24) * x28); }
template<> inline double xi<17>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; double x216 = x28 * x28; return x * (9.9999999999991586e-01 + -3.3333333331539649e-01 * x2 + (1.9999999885901307e-01 + -1.4285710866597597e-01 * x2) * x22 + (1.1111052308663949e-01 + -9.0902614649300906e-02 * x2 + (7.6874146803011537e-02 + -6.6400996804419749e-02 * x2) * x22) * x24 + (5.7751981894433748e-02 + -4.9340121701547388e-02 * x2 + (3.9762413378979828e-02 + -2.8629853703447561e-02 * x2) * x22 + (1.7310365496208254e-02 + -8.2096187905995513e-03 * x2 + (2.8093457170694520e-03 + -6.0940031581488599e-04 * x2) * x22) * x24) * x28 + 6.2436108681934464e-05 * x216); }
template<> inline double xi<18>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; double x216 = x28 * x28; return x * (9.9999999999998558e-01 + -3.3333333332990472e-01 * x2 + (1.9999999975665970e-01 + -1.4285713471275124e-01 * x2) * x22 + (1.1111095442381502e-01 + -9.0907156394328137e-02 * x2 + (7.6906649156526998e-02 + -6.6566100086412207e-02 * x2) * x22) * x24 + (5.8364603332581473e-02 + -5.1030977240766378e-02 * x2 + (4.3267805093829844e-02 + -3.4100063401437874e-02 * x2) * x22 + (2.3697852281270090e-02 + -1.3702546920513891e-02 * x2 + (6.1837955998122361e-03 + -2.0101395373502504e-03 * x2) * x22) * x24) * x28 + (4.1436296328711983e-04 + -4.0407586855442829e-05 * x2) * x216); }
template<> inline double xi<19>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; double x216 = x28 * x28; return x * (9.9999999999999753e-01 + -3.3333333333268159e-01 * x2 + (1.9999999994868435e-01 + -1.4285714095010484e-01 * x2) * x22 + (1.1111107031789423e-01 + -9.0908529899004480e-02 * x2 + (7.6917758114345896e-02 + -6.6630215499641454e-02 * x2) * x22) * x24 + (5.8636732908059056e-02 + -5.1897548887283423e-02 * x2 + (4.5363655637590359e-02 + -3.7971289915217645e-02 * x2) * x22 + (2.9153440169624127e-02 + -1.9520100422746680e-02 * x2 + (1.0795336607971266e-02 + -4.6423188222075775e-03 * x2) * x22) * x24) * x28 + (1.4364372066614275e-03 + -2.8202613038791620e-04 * x2 + 2.6236345895738590e-05 * x22) * x216); }
template<> inline double xi<20>(double x) { double x2 = x * x; double x22 = x2 * x2; double x24 = x22 * x22; double x28 = x24 * x24; double x216 = x28 * x28; return x * (9.9999999999999958e-01 + -3.3333333333321006e-01 * x2 + (1.9999999998928794e-01 + -1.4285714241743987e-01 * x2) * x22 + (1.1111110071008497e-01 + -9.0908932480142952e-02 * x2 + (7.6921409987427713e-02 + -6.6653958158564996e-02 * x2) * x22) * x24 + (5.8750881194559933e-02 + -5.2312166516808036e-02 * x2 + (4.6517575683483294e-02 + -4.0451703413633348e-02 * x2) * x22 + (3.3281457966199757e-02 + -2.4821296843759711e-02 * x2 + (1.5995635936942689e-02 + -8.4651086471683459e-03 * x2) * x22) * x24) * x28 + (3.4738841849962415e-03 + -1.0251712186904009e-03 * x2 + (1.9211551727884094e-04 + -1.7084743394938366e-05 * x2) * x22) * x216); }

}  // namespace ecef2geo
}  // namespace priv

#endif // ECEF2GEO_XI_HPP
