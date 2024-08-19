/*
 * Generated. DO NOT EDIT!
 */

#ifndef ECEF2GEO_OMEGA_HPP
#define ECEF2GEO_OMEGA_HPP

namespace ecef2geo {

// Minimax approximation ranges.
constexpr double OMEGA_H_MIN = -5.00000000000000000e+03;
constexpr double OMEGA_H_MAX = 1.00000000000000000e+05;
constexpr double OMEGA_H_0 = 0.00000000000000000e+00;
// For ensuring compile-time removal of h_0.
#define H_0_IS_ZERO

namespace priv {

/** Polynomial additive latitude correction
 *
 * Nilsson, J.-O. Minimax polynomial ECEF to geodetic coordinate transformation
 * approximations. doi: https://doi.org/10.48550/arXiv.2405.06572
 *
 * @paramt N Number of Fourier coefficients in the approximation.
 * @paramt M Number of coefficient for the altitude dependency.
 * @param u h_c.
 * @param v t^2.
 * @return Minimax approximation of omega.
 */
template<int N,int M> inline double omega(double u, double v) = delete;  // Only allow provided specializations.
template<> inline double omega<0,0>(double u, double v) { return 0.0; }
template<> inline double omega<0,1>(double u, double v) { return 0.0; }
template<> inline double omega<0,2>(double u, double v) { return 0.0; }
template<> inline double omega<0,3>(double u, double v) { return 0.0; }
template<> inline double omega<0,4>(double u, double v) { return 0.0; }
template<> inline double omega<0,5>(double u, double v) { return 0.0; }
template<> inline double omega<0,6>(double u, double v) { return 0.0; }
template<> inline double omega<1,0>(double u, double v) { return 0.0066677813753770136; }
template<> inline double omega<1,1>(double u, double v) { return 0.013335202317269038 + -1.0394079474920395e-09 * u; }
template<> inline double omega<1,2>(double u, double v) { double u2 = u * u; return 0.020002558556087308 + -3.1182288863767869e-09 * u + 1.6203149568692841e-16 * u2; }
template<> inline double omega<1,3>(double u, double v) { double u2 = u * u; return 0.026669822295348564 + -6.2364044268901539e-09 * u + (6.4812710099485553e-16 + -2.5258702517139241e-23 * u) * u2; }
template<> inline double omega<1,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; return 0.033336965747610005 + -1.0393888399376689e-08 * u + (1.6203047026245014e-15 + -1.2629352805875582e-22 * u) * u2 + 3.9375023366464455e-30 * u4; }
template<> inline double omega<1,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; return 0.040003961147859715 + -1.5590621657559157e-08 * u + (3.2405673447681642e-15 + -3.7887712753666839e-22 * u) * u2 + (2.3624958916061237e-29 + -6.1380271062803481e-37 * u) * u4; }
template<> inline double omega<1,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; return 0.046670780678684999 + -2.1826532023634562e-08 * u + (5.6708991630794645e-15 + -8.8403345218192844e-22 * u) * u2 + (8.2686439383985845e-29 + -4.2965962847706622e-36 * u + 9.568304070813024e-44 * u2) * u4; }
template<> inline double omega<2,0>(double u, double v) { return 0.0067010483124497356 + -6.6533874145443769e-05 * v; }
template<> inline double omega<2,1>(double u, double v) { return 0.013446184736230014 + -1.0515236264437181e-09 * u + (-0.00022196483792195034 + 2.4231357903357331e-11 * u) * v; }
template<> inline double omega<2,2>(double u, double v) { double u2 = u * u; return 0.020235710992772528 + -3.1684362320910185e-09 * u + 1.6500052340136923e-16 * u2 + (-0.00046630487337043957 + 1.0041469142846286e-10 * u + -5.9380554288815973e-18 * u2) * v; }
template<> inline double omega<2,3>(double u, double v) { double u2 = u * u; return 0.027069596697081113 + -6.3645384419626669e-09 * u + (6.6324422388573312e-16 + -2.5889948265686474e-23 * u) * u2 + (-0.00079954880346509607 + 2.5626803014502601e-10 * u + (-3.0234245781755037e-17 + 1.2624914970944679e-24 * u) * u2) * v; }
template<> inline double omega<2,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; return 0.033947812057822907 + -1.0653642380896784e-08 * u + (1.6661996102540888e-15 + -1.3012339932444736e-22 * u) * u2 + 4.0621584784028504e-30 * u4 + (-0.0012216926204258079 + 5.1950796304019089e-10 * u + (-9.1789815259174709e-17 + 7.6597425313830808e-24 * u) * u2 + -2.4931228351280943e-31 * u4) * v; }
template<> inline double omega<2,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; return 0.040870326813284015 + -1.6049545875811873e-08 * u + (3.3485605650039408e-15 + -3.9238754574022328e-22 * u) * u2 + (2.4504159330539716e-29 + -6.3732742083748206e-37 * u) * u4 + (-0.0017327313308486047 + 9.1784843650543299e-10 * u + (-2.1598644047155308e-16 + 2.7020836407109731e-23 * u) * u2 + (-1.7584008289569594e-30 + 4.7049420418894447e-38 * u) * u4) * v; }
template<> inline double omega<2,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; return 0.047837110152839389 + -2.2566032090161163e-08 * u + (5.8882417453288786e-15 + -9.2027268235267315e-22 * u) * u2 + (8.6223031229448686e-29 + -4.4858234100309645e-36 * u + 9.9988174869285852e-44 * u2) * u4 + (-0.0023326589483087775 + 1.4790001330532036e-09 * u + (-4.3468516449882893e-16 + 7.2478460341489239e-23 * u) * u2 + (-7.07318369092568e-30 + 3.7845425052060462e-37 * u + -8.6102683223112143e-45 * u2) * u4) * v; }
template<> inline double omega<3,0>(double u, double v) { double v2 = v * v; return 0.0067012258007855211 + -6.7480478602966758e-05 * v + 9.4660445752299873e-07 * v2; }
template<> inline double omega<3,1>(double u, double v) { double v2 = v * v; return 0.013447020085588264 + -1.0516261888232214e-09 * u + (-0.00022642003449928512 + 2.4778357260709133e-11 * u) * v + (4.4551965773347914e-06 + -5.4699935735179955e-13 * u) * v2; }
template<> inline double omega<3,2>(double u, double v) { double u2 = u * u; double v2 = v * v; return 0.02023800884030541 + -3.1689947923483929e-09 * u + 1.6503606600562022e-16 * u2 + (-0.00047856006021247579 + 1.0339367946779417e-10 * u + -6.1276159848868637e-18 * u2) * v + (1.2255186842036269e-05 + -2.9789880393313049e-12 * u + 1.8956055600526669e-19 * u2) * v2; }
template<> inline double omega<3,3>(double u, double v) { double u2 = u * u; double v2 = v * v; return 0.027074485817535208 + -6.3663089057532568e-09 * u + (6.6346869224010302e-16 + -2.5899765348800763e-23 * u) * u2 + (-0.00082562411255362067 + 2.6571050369483905e-10 * u + (-3.1431410338394555e-17 + 1.3148492737039995e-24 * u) * u2) * v + (2.60753090885247e-05 + -9.4424735498130361e-12 * u + (1.1971645566395157e-18 + -5.2357776609531562e-26 * u) * u2) * v2; }
template<> inline double omega<3,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; return 0.033956745360127084 + -1.0657934714234198e-08 * u + (1.6670137890290881e-15 + -1.3019450299569873e-22 * u) * u2 + 4.0645469293734972e-30 * u4 + (-0.0012693368993814277 + 5.4240040750639846e-10 * u + (-9.6132102059170536e-17 + 8.038962111390361e-24 * u) * u2 + -2.6205068868959098e-31 * u4) * v + (4.7644278955619766e-05 + -2.2892444466207578e-11 * u + (4.3422867999958359e-18 + -3.7921958000728011e-25 * u) * u2 + 1.2738405176781577e-32 * u4) * v2; }
template<> inline double omega<3,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; return 0.040885081311572244 + -1.6058375673770013e-08 * u + (3.3507894589328667e-15 + -3.9267919049065346e-22 * u) * u2 + (2.4523737735781125e-29 + -6.3786335898484567e-37 * u) * u4 + (-0.0018114219883858234 + 9.6494069228217699e-10 * u + (-2.2787387475915666e-16 + 2.8576275076070839e-23 * u) * u2 + (-1.862818990244474e-30 + 4.9907757204833766e-38 * u) * u4) * v + (7.8690657537218718e-05 + -4.7092255776744055e-11 * u + (1.1887434287603587e-17 + -1.5554386689611109e-24 * u) * u2 + (1.0441816128751447e-31 + -2.8583367859393255e-39 * u) * u4) * v2; }
template<> inline double omega<3,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; return 0.047859786931290428 + -2.2582272121727575e-08 * u + (5.8933586456203936e-15 + -9.2116461405762506e-22 * u) * u2 + (8.6312793612201355e-29 + -4.4907356932902443e-36 * u + 1.0010187728496843e-43 * u2) * u4 + (-0.0024536017667143336 + 1.5656136347407201e-09 * u + (-4.6197529938690531e-16 + 7.7235429434566839e-23 * u) * u2 + (-7.5519163989399298e-30 + 4.0465309457009936e-37 * u + -9.2166812059515464e-45 * u2) * u4) * v + (0.00012094281840555624 + -8.6613501687516475e-11 * u + (2.7290134888076368e-17 + -4.7569690930775963e-24 * u) * u2 + (4.7873270801424953e-31 + -2.6198844049494753e-38 * u + 6.0641288364033224e-46 * u2) * u4) * v2; }
template<> inline double omega<4,0>(double u, double v) { double v2 = v * v; return 0.0067012267744125092 + -6.7490214872849025e-05 * v + (9.6997150524041982e-07 + -1.5578031811614025e-08 * v) * v2; }
template<> inline double omega<4,1>(double u, double v) { double v2 = v * v; return 0.013447026025634832 + -1.0516269631274677e-09 * u + (-0.00022647943496496668 + 2.4786100303171601e-11 * u) * v + (4.5977576949705561e-06 + -5.6558265926172858e-13 * u + (-9.5040745090509506e-08 + 1.2388867939952638e-14 * u) * v) * v2; }
template<> inline double omega<4,2>(double u, double v) { double u2 = u * u; double v2 = v * v; return 0.020238029191756247 + -3.1690000600628975e-09 * u + 1.6503641624484183e-16 * u2 + (-0.00047876357472084296 + 1.0344635661283742e-10 * u + -6.1311183771030082e-18 * u2) * v + (1.2743621662117424e-05 + -3.1054131874350789e-12 * u + 1.9796629732401303e-19 * u2 + (-3.256232133874365e-07 + 8.4283432069182681e-14 * u + -5.6038275458308947e-21 * u2) * v) * v2; }
template<> inline double omega<4,3>(double u, double v) { double u2 = u * u; double v2 = v * v; return 0.027074537946419158 + -6.3663290353839142e-09 * u + (6.6347135933514977e-16 + -2.5899885739079588e-23 * u) * u2 + (-0.00082614540139308665 + 2.6591180000141378e-10 * u + (-3.1458081288862306e-17 + 1.3160531764922644e-24 * u) * u2) * v + (2.7326402303242933e-05 + -9.9255846855924077e-12 * u + (1.2611748377621174e-18 + -5.524714330136761e-26 * u) * u2 + (-8.3406214314548869e-07 + 3.2207409051958092e-13 * u + (-4.2673520748401045e-20 + 1.9262444612240382e-27 * u) * u2) * v) * v2; }
template<> inline double omega<4,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; return 0.033956857023526403 + -1.0657991968432448e-08 * u + (1.6670251372934867e-15 + -1.3019552559156135e-22 * u) * u2 + 4.0645820900349371e-30 * u4 + (-0.0012704535333745719 + 5.4297294948889652e-10 * u + (-9.6245584703157008e-17 + 8.0491880700165421e-24 * u) * u2 + -2.6240229530399466e-31 * u4) * v + (5.0324200539165746e-05 + -2.4266545224202978e-11 * u + (4.6146451455633588e-18 + -4.0376188071011581e-25 * u) * u2 + 1.3582261051350419e-32 * u4 + (-1.7866143890306541e-06 + 9.1606717199693324e-13 * u + (-1.8157223037834837e-19 + 1.6361533801890488e-26 * u) * u2 + -5.6257058304589458e-34 * u4) * v) * v2; }
template<> inline double omega<4,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; return 0.040885293127108015 + -1.6058510993921203e-08 * u + (3.3508251470438082e-15 + -3.9268400744834055e-22 * u) * u2 + (2.4524068646191565e-29 + -6.3787257970306808e-37 * u) * u4 + (-0.0018135401437434969 + 9.6629389379406606e-10 * u + (-2.2823075586857513e-16 + 2.8624444652941675e-23 * u) * u2 + (-1.8661280943488577e-30 + 4.999996438705822e-38 * u) * u4) * v + (8.3774230395635378e-05 + -5.0339939405277698e-11 * u + (1.2743948950207908e-17 + -1.671045653451113e-24 * u) * u2 + (1.123600111380355e-31 + -3.0796340232780235e-39 * u) * u4 + (-3.3890485722777727e-06 + 2.1651224190224279e-12 * u + (-5.7100977506954692e-19 + 7.7071322993334787e-26 * u) * u2 + (-5.294566567014007e-33 + 1.4753149155913206e-40 * u) * u4) * v) * v2; }
template<> inline double omega<4,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; return 0.047860154846244451 + -2.2582553452106821e-08 * u + (5.8934512386514108e-15 + -9.2118125900315352e-22 * u) * u2 + (8.6314507420881839e-29 + -4.4908311451451016e-36 * u + 1.0010411766932686e-43 * u2) * u4 + (-0.0024572809162545423 + 1.5684269385331804e-09 * u + (-4.6290122969708122e-16 + 7.7401878889851026e-23 * u) * u2 + (-7.5690544857447023e-30 + 4.0560761311867104e-37 * u + -9.2390850495359153e-45 * u2) * u4) * v + (0.0001297727773020576 + -9.3365430789421105e-11 * u + (2.9512367632498714e-17 + -5.1564477857596539e-24 * u) * u2 + (5.1986411634570113e-31 + -2.8489688566066765e-38 * u + 6.6018210824281669e-46 * u2) * u4 + (-5.8866392643342399e-06 + 4.5012860679364211e-12 * u + (-1.4814884962815641e-18 + 2.66319128454705e-25 * u) * u2 + (-2.7420938887634368e-32 + 1.5272296777146727e-39 * u + -3.5846149734989648e-47 * u2) * u4) * v) * v2; }
template<> inline double omega<5,0>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.0067012267798402889 + -6.749030171732508e-05 * v + (9.7033625203985761e-07 + -1.6133836458376546e-08 * v) * v2 + 2.7790232338126106e-10 * v4; }
template<> inline double omega<5,1>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.013447026066401196 + -1.0516269686372191e-09 * u + (-0.00022648008722678744 + 2.4786188459191202e-11 * u) * v + (4.600497194617681e-06 + -5.6595291454405164e-13 * u + (-9.9215220743271937e-08 + 1.2953066465397291e-14 * u) * v) * v2 + (2.0872378263812176e-09 + -2.8209926272232652e-16 * u) * v4; }
template<> inline double omega<5,2>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.020238029359570701 + -3.1690001051858808e-09 * u + 1.6503641933251515e-16 * u2 + (-0.00047876625975213391 + 1.0344707858057174e-10 * u + -6.1311677798763597e-18 * u2) * v + (1.2754898793539329e-05 + -3.1084454519192448e-12 * u + 1.9817378897208824e-19 * u2 + (-3.428074136493872e-07 + 8.8904025568864013e-14 * u + -5.9200052952788245e-21 * u2) * v) * v2 + (8.5921001309753337e-09 + -2.3102967498406646e-15 * u + 1.5808887472396487e-22 * u2) * v4; }
template<> inline double omega<5,3>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.02707453845132074 + -6.3663292381588322e-09 * u + (6.6347138699958437e-16 + -2.5899887016159891e-23 * u) * u2 + (-0.0008261534798184089 + 2.6591504440010536e-10 * u + (-3.1458523919816397e-17 + 1.3160736097771083e-24 * u) * u2) * v + (2.7360331689596588e-05 + -9.939211160097003e-12 * u + (1.2630338877693007e-18 + -5.5332963097711163e-26 * u) * u2 + (-8.8576406520820212e-07 + 3.4283824214563153e-13 * u + (-4.5506358854585118e-20 + 2.0570174842237393e-27 * u) * u2) * v) * v2 + (2.5850961031356672e-08 + -1.0382075813025301e-14 * u + (1.4164190530920376e-21 + -6.5386511499850515e-29 * u) * u2) * v4; }
template<> inline double omega<5,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.033956858269743512 + -1.0657992633478044e-08 * u + (1.6670252730549704e-15 + -1.3019553810283226e-22 * u) * u2 + 4.0645825278533089e-30 * u4 + (-0.0012704734728483523 + 5.4298359021842561e-10 * u + (-9.6247756886895796e-17 + 8.0493882503509578e-24 * u) * u2 + -2.6240930039793388e-31 * u4) * v + (5.0407946329043319e-05 + -2.4311236288225156e-11 * u + (4.6237683172662731e-18 + -4.046026381146591e-25 * u) * u2 + 1.3611682445895033e-32 * u4 + (-1.9142270212250502e-06 + 9.8416784098310913e-13 * u + (-1.9547420630659842e-19 + 1.7642687942146882e-26 * u) * u2 + -6.074031842567329e-34 * u4) * v) * v2 + (6.3806316097197935e-08 + -3.4050334493087929e-14 * u + (6.9509879641250328e-21 + -6.4057707012819846e-28 * u) * u2 + 2.2416300605419168e-35 * u4) * v4; }
template<> inline double omega<5,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.040885295808229409 + -1.6058512777438019e-08 * u + (3.3508256315289557e-15 + -3.9268407432252944e-22 * u) * u2 + (2.4524073321306264e-29 + -6.3787271181122606e-37 * u) * u4 + (-0.0018135830416858168 + 9.6632243006313259e-10 * u + (-2.2823850763093337e-16 + 2.8625514639963372e-23 * u) * u2 + (-1.8662028961840396e-30 + 5.0002078117585103e-38 * u) * u4) * v + (8.3954401753378994e-05 + -5.0459791735357369e-11 * u + (1.2776506352112479e-17 + -1.6755395989422308e-24 * u) * u2 + (1.1267417884579888e-31 + -3.0885116914909313e-39 * u) * u4 + (-3.6635954031251799e-06 + 2.3477545410485884e-12 * u + (-6.2062105416222606e-19 + 8.391923993218097e-26 * u) * u2 + (-5.7732983121772601e-33 + 1.6105936693118173e-40 * u) * u4) * v) * v2 + (1.3727341542370365e-07 + -9.1316061013080278e-14 * u + (2.4805639546339552e-20 + -3.4239584694230862e-27 * u) * u2 + (2.3936587258162673e-34 + -6.7639376860248341e-42 * u) * u4) * v4; }
template<> inline double omega<5,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.047860160059700063 + -2.2582557604286886e-08 * u + (5.893452646282074e-15 + -9.2118151775814239e-22 * u) * u2 + (8.6314534530227597e-29 + -4.4908326761490727e-36 * u + 1.0010415401427854e-43 * u2) * u4 + (-0.0024573643315443068 + 1.5684933734142425e-09 * u + (-4.6292375178769741e-16 + 7.7406018969673191e-23 * u) * u2 + (-7.5694882352769732e-30 + 4.0563210918220539e-37 * u + -9.2396665687627159e-45 * u2) * u4) * v + (0.00013012312151906917 + -9.3644457289882182e-11 * u + (2.9606960413086595e-17 + -5.1738361210127322e-24 * u) * u2 + (5.216858643812416e-31 + -2.8592572032910933e-38 * u + 6.6262448899538062e-46 * u2) * u4 + (-6.420497118828065e-06 + 4.926469306734239e-12 * u + (-1.6256298762250038e-18 + 2.928156393165386e-25 * u) * u2 + (-3.0196935894172255e-32 + 1.6840044843343569e-39 * u + -3.956787278651553e-47 * u2) * u4) * v) * v2 + (2.6692892724691259e-07 + -2.1259161939890916e-13 * u + (7.2070689971719853e-20 + -1.324825543091681e-26 * u) * u2 + (1.3879985032689433e-33 + -7.8387403309842138e-41 * u + 1.8608615257629417e-48 * u2) * u4) * v4; }
template<> inline double omega<6,0>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.0067012267798708756 + -6.7490302431005646e-05 * v + (9.7034081959543981e-07 + -1.6145581601301731e-08 * v) * v2 + (2.9095248218702237e-10 + -5.2200635223045176e-12 * v) * v4; }
template<> inline double omega<6,1>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.013447026066674334 + -1.0516269686750374e-09 * u + (-0.000226480093600012 + 2.4786189341618833e-11 * u) * v + (4.6005379832549828e-06 + -5.6595856208087875e-13 * u + (-9.9320105810619268e-08 + 1.296758870295277e-14 * u) * v) * v2 + (2.2037767901004768e-09 + -2.982350822284147e-16 * u + (-4.6615585487703673e-11 + 6.4543278024352727e-18 * u) * v) * v4; }
template<> inline double omega<6,2>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.020238029360886305 + -3.1690001055487389e-09 * u + 1.6503641935785068e-16 * u2 + (-0.00047876629044955329 + 1.0344708704725596e-10 * u + -6.1311683710383436e-18 * u2) * v + (1.2755095257023495e-05 + -3.1084996386982771e-12 * u + 1.9817757240878402e-19 * u2 + (-3.4331260546581518e-07 + 8.9043363000661733e-14 * u + -5.9297341324965561e-21 * u2) * v) * v2 + (9.1534243714509172e-09 + -2.4651161185048045e-15 * u + 1.6889869385477765e-22 * u2 + (-2.2452969619023332e-10 + 6.1927747465655956e-17 * u + -4.323927652325117e-24 * u2) * v) * v4; }
template<> inline double omega<6,3>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.027074538455879524 + -6.3663292400384943e-09 * u + (6.634713872613789e-16 + -2.5899887028447025e-23 * u) * u2 + (-0.00082615358619002321 + 2.6591508825887243e-10 * u + (-3.1458530028355177e-17 + 1.316073896476898e-24 * u) * u2) * v + (2.7361012467927818e-05 + -9.9394918562063401e-12 * u + (1.2630729824174771e-18 + -5.5334797976365962e-26 * u) * u2 + (-8.8751463805993047e-07 + 3.4356003214106887e-13 * u + (-4.5606887949896142e-20 + 2.0617357436217837e-27 * u) * u2) * v) * v2 + (2.7796041977721567e-08 + -1.1184064696844569e-14 * u + (1.5281180478820613e-21 + -7.062902194212217e-29 * u) * u2 + (-7.7803237854595732e-10 + 3.2079555352770725e-16 * u + (-4.4679597916009476e-23 + 2.0970041769086604e-30 * u) * u2) * v) * v4; }
template<> inline double omega<6,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.033956858282520062 + -1.0657992640482159e-08 * u + (1.6670252745150658e-15 + -1.3019553823965527e-22 * u) * u2 + 4.0645825327067229e-30 * u4 + (-0.0012704737709678274 + 5.4298375364780816e-10 * u + (-9.6247790955787215e-17 + 8.0493914428879671e-24 * u) * u2 + -2.6240941364427088e-31 * u4) * v + (5.0409854293683779e-05 + -2.4312282236272959e-11 * u + (4.6239863581713787e-18 + -4.0462307035152413e-25 * u) * u2 + 1.3612407222451911e-32 * u4 + (-1.9191332160148131e-06 + 9.8685742167745805e-13 * u + (-1.9603488291972698e-19 + 1.7695227979799916e-26 * u) * u2 + -6.0926689540299309e-34 * u4) * v) * v2 + (6.9257643641378892e-08 + -3.703875748680897e-14 * u + (7.5739619787123119e-21 + -6.9895488974267973e-28 * u) * u2 + 2.4487090767930503e-35 * u4 + (-2.1805310176723792e-09 + 1.1953691974884166e-15 * u + (-2.4918960583491195e-22 + 2.335112784579253e-29 * u) * u2 + -8.2831606500453459e-37 * u4) * v) * v4; }
template<> inline double omega<6,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.040885295839059019 + -1.6058512798514048e-08 * u + (3.3508256373764818e-15 + -3.9268407514331456e-22 * u) * u2 + (2.452407337947089e-29 + -6.3787271347333624e-37 * u) * u4 + (-0.0018135837610434001 + 9.6632292183719524e-10 * u + (-2.2823864407320369e-16 + 2.8625533791616239e-23 * u) * u2 + (-1.8662042533586911e-30 + 5.0002116900156677e-38 * u) * u4) * v + (8.3959005641911052e-05 + -5.0462939089357817e-11 * u + (1.2777379582642763e-17 + -1.6756621695205854e-24 * u) * u2 + (1.1268286476356885e-31 + -3.0887598999489422e-39 * u) * u4 + (-3.6754339736361957e-06 + 2.3558477370497502e-12 * u + (-6.2286650409724249e-19 + 8.4234421419378366e-26 * u) * u2 + (-5.7956335293000243e-33 + 1.6169761725178128e-40 * u) * u4) * v) * v2 + (1.5042738265816568e-07 + -1.0030850101437109e-13 * u + (2.7300583918580043e-20 + -3.7741601218646458e-27 * u) * u2 + (2.6418278049580857e-34 + -7.473104708913223e-42 * u) * u4 + (-5.2615868937848073e-09 + 3.5969760005163239e-15 * u + (-9.9797774889619664e-22 + 1.4008066097662379e-28 * u) * u2 + (-9.926763165672736e-36 + 2.8366680915535596e-43 * u) * u4) * v) * v4; }
template<> inline double omega<6,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.047860160126233016 + -2.2582557658758688e-08 * u + (5.8934526651450343e-15 + -9.2118152128425794e-22 * u) * u2 + (8.6314534904693032e-29 + -4.4908326975342654e-36 * u + 1.0010415452670768e-43 * u2) * u4 + (-0.0024573658839798789 + 1.5684946444229856e-09 * u + (-4.6292419192342094e-16 + 7.7406101245702341e-23 * u) * u2 + (-7.5694969728037387e-30 + 4.0563260817003679e-37 * u + -9.2396785254430326e-45 * u2) * u4) * v + (0.00013013305710672752 + -9.3652591745837749e-11 * u + (2.9609777281717268e-17 + -5.1743626875992793e-24 * u) * u2 + (5.2174178455253861e-31 + -2.8595765555032261e-38 * u + 6.6270101174940671e-46 * u2) * u4 + (-6.4460457728067195e-06 + 4.9473864791914186e-12 * u + (-1.6328732527038827e-18 + 2.9416966768194697e-25 * u) * u2 + (-3.0340730620364553e-32 + 1.6922163983606351e-39 * u + -3.9764645582582794e-47 * u2) * u4) * v) * v2 + (2.9531632055652862e-07 + -2.3583292212910806e-13 * u + (8.0118886059363075e-20 + -1.4752731392481644e-26 * u) * u2 + (1.547770421260386e-33 + -8.7511752227929099e-41 * u + 2.0794979658376787e-48 * u2) * u4 + (-1.1354957323846409e-08 + 9.2965210920795608e-15 * u + (-3.21927843505729e-21 + 6.0179038462593335e-28 * u) * u2 + (-6.3908767196577137e-35 + 3.6497395672347853e-42 * u + -8.7454576029894789e-50 * u2) * u4) * v) * v4; }
template<> inline double omega<7,0>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.006701226779871049 + -6.7490302436564202e-05 * v + (9.7034086962253382e-07 + -1.6145772180707711e-08 * v) * v2 + (2.9130187776465272e-10 + -5.524990571872839e-12 * v + 1.0164234985610733e-13 * v2) * v4; }
template<> inline double omega<7,1>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.013447026066676132 + -1.0516269686752907e-09 * u + (-0.00022648009365757292 + 2.4786189349727239e-11 * u) * v + (4.600538501303193e-06 + -5.6595863505654166e-13 * u + (-9.9322079327608881e-08 + 1.2967866705478107e-14 * u) * v) * v2 + (2.2073949045814152e-09 + -2.9874475352486707e-16 * u + (-4.9773212671068153e-11 + 6.8991318429755161e-18 * u) * v + (1.0525423944548254e-12 + -1.4826801351341438e-19 * u) * v2) * v4; }
template<> inline double omega<7,2>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.020238029360896241 + -3.1690001055515294e-09 * u + 1.6503641935804843e-16 * u2 + (-0.00047876629076749444 + 1.034470871365511e-10 * u + -6.1311683773665542e-18 * u2) * v + (1.2755098118494028e-05 + -3.1085004423545023e-12 * u + 1.9817762936268572e-19 * u2 + (-3.4332350630593771e-07 + 8.9046424548186028e-14 * u + -5.9299510997410841e-21 * u2) * v) * v2 + (9.1734092450088263e-09 + -2.4707289556326779e-15 * u + 1.6929646713641272e-22 * u2 + (-2.4197104038622714e-10 + 6.6826223504527349e-17 * u + -4.6710752435702608e-24 * u2) * v + (5.8137813986646008e-12 + -1.6328253462904638e-18 * u + 1.1571586374838119e-25 * u2) * v2) * v4; }
template<> inline double omega<7,3>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.027074538455918552 + -6.3663292400548915e-09 * u + (6.6347138726369786e-16 + -2.5899887028557251e-23 * u) * u2 + (-0.00082615358743896319 + 2.6591508878358784e-10 * u + (-3.1458530102562534e-17 + 1.3160739000041161e-24 * u) * u2) * v + (2.7361023708388118e-05 + -9.9394965786451695e-12 * u + (1.2630736502837104e-18 + -5.5334829721328621e-26 * u) * u2 + (-8.8755745886108134e-07 + 3.435780223842309e-13 * u + (-4.5609432202213442e-20 + 2.0618566768128568e-27 * u) * u2) * v) * v2 + (2.7874546779831435e-08 + -1.1217046809308308e-14 * u + (1.5327825104637807e-21 + -7.085073279242314e-29 * u) * u2 + (-8.4654566038729769e-10 + 3.4957994258696981e-16 * u + (-4.8750401623691814e-23 + 2.2904972826258691e-30 * u) * u2) * v + (2.2837760613780126e-11 + -9.594796353087531e-18 * u + (1.3569345692274461e-24 + -6.4497701905736243e-32 * u) * u2) * v2) * v4; }
template<> inline double omega<7,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.033956858282642735 + -1.0657992640550716e-08 * u + (1.6670252745295816e-15 + -1.3019553824103309e-22 * u) * u2 + 4.0645825327561236e-30 * u4 + (-0.0012704737748933813 + 5.4298375584161584e-10 * u + (-9.6247791420295202e-17 + 8.0493914869780477e-24 * u) * u2 + -2.6240941522508999e-31 * u4) * v + (5.0409889623668919e-05 + -2.4312301980542879e-11 * u + (4.6239905387432534e-18 + -4.0462346716224499e-25 * u) * u2 + 1.3612421449823831e-32 * u4 + (-1.9192678064344007e-06 + 9.8693263794381655e-13 * u + (-1.9605080890782172e-19 + 1.7696739639688808e-26 * u) * u2 + -6.0932109491506457e-34 * u4) * v) * v2 + (6.9504392743956241e-08 + -3.7176653975132853e-14 * u + (7.6031596235526873e-21 + -7.0172626620564421e-28 * u) * u2 + 2.4586456540061538e-35 * u4 + (-2.3958756890126109e-09 + 1.3157152236619882e-15 * u + (-2.746711867865129e-22 + 2.5769783668016065e-29 * u) * u2 + -9.1503528431888913e-37 * u4) * v + (7.1781557113410595e-11 + -4.011534205785723e-17 * u + (8.4938603172003096e-24 + -8.0621860740784527e-31 * u) * u2 + 2.8906406438118198e-38 * u4) * v2) * v4; }
template<> inline double omega<7,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.04088529583938779 + -1.6058512798743253e-08 * u + (3.3508256374410852e-15 + -3.9268407515250062e-22 * u) * u2 + (2.4524073380128906e-29 + -6.3787271349231128e-37 * u) * u4 + (-0.0018135837715640589 + 9.6632292917172942e-10 * u + (-2.2823864614052162e-16 + 2.8625534085570648e-23 * u) * u2 + (-1.8662042744151284e-30 + 5.0002117507358349e-38 * u) * u4) * v + (8.3959100327840857e-05 + -5.0463005100166592e-11 * u + (1.2777398188503717e-17 + -1.6756648151102934e-24 * u) * u2 + (1.1268305427150402e-31 + -3.0887653647639996e-39 * u) * u4 + (-3.6757946819402147e-06 + 2.3560992067974613e-12 * u + (-6.2293738356754512e-19 + 8.4244499856361661e-26 * u) * u2 + (-5.7963554642911529e-33 + 1.6171843559485748e-40 * u) * u4) * v) * v2 + (1.510886812155338e-07 + -1.0076952888517477e-13 * u + (2.7430529614134862e-20 + -3.7926372563340241e-27 * u) * u2 + (2.6550632797954522e-34 + -7.5112716712196043e-42 * u) * u4 + (-5.8387201802151793e-09 + 3.9993275968540849e-15 * u + (-1.111384901380404e-21 + 1.5620616014989919e-28 * u) * u2 + (-1.1081859151479268e-35 + 3.169761580772885e-43 * u) * u4) * v + (1.9237776214345757e-10 + -1.3411719877925391e-16 * u + (3.7802384161402447e-23 + -5.3751663910918058e-30 * u) * u2 + (3.8503199526884365e-37 + -1.1103116307310838e-44 * u) * u4) * v2) * v4; }
template<> inline double omega<7,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.047860160127014342 + -2.2582557659411201e-08 * u + (5.8934526653746141e-15 + -9.2118152132773543e-22 * u) * u2 + (8.6314534909360319e-29 + -4.4908326978032404e-36 * u + 1.00104154533203e-43 * u2) * u4 + (-0.0024573659089823312 + 1.5684946653033398e-09 * u + (-4.6292419926999548e-16 + 7.7406102636980474e-23 * u) * u2 + (-7.5694971221567537e-30 + 4.0563261677723679e-37 * u + -9.2396787332926614e-45 * u2) * u4) * v + (0.00013013328212880239 + -9.3652779669025828e-11 * u + (2.9609843400888225e-17 + -5.1743752091025225e-24 * u) * u2 + (5.2174312872967964e-31 + -2.8595843019831626e-38 * u + 6.6270288239607048e-46 * u2) * u4 + (-6.4469029997585337e-06 + 4.9481023770507638e-12 * u + (-1.6331251352598994e-18 + 2.9421736864668261e-25 * u) * u2 + (-3.0345851295187451e-32 + 1.692511502358207e-39 * u + -3.9771771855587394e-47 * u2) * u4) * v) * v2 + (2.968879033015212e-07 + -2.3714540153790882e-13 * u + (8.0580670745393508e-20 + -1.4840183161163701e-26 * u) * u2 + (1.5571583251023658e-33 + -8.8052776223477564e-41 * u + 2.0925627996794424e-48 * u2) * u4 + (-1.2726520446749035e-08 + 1.0441957667032937e-14 * u + (-3.6222905246838521e-21 + 6.7811192820300256e-28 * u) * u2 + (-7.2101846913213811e-35 + 4.1219059633498055e-42 * u + -9.8856612837252178e-50 * u2) * u4) * v + (4.5718770763420863e-10 + -3.818121916511254e-16 * u + (1.3433736320885411e-22 + -2.5440514525689721e-29 * u) * u2 + (2.7310265722122229e-36 + -1.5738879870500664e-43 * u + 3.8006789357857993e-51 * u2) * u4) * v2) * v4; }
template<> inline double omega<8,0>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.0067012267798710499 + -6.7490302436605889e-05 * v + (9.7034087012269428e-07 + -1.6145774800595489e-08 * v) * v2 + (2.9130886413205713e-10 + -5.5348985111008952e-12 * v + (1.0875574212240394e-13 + -2.0323977903704622e-15 * v) * v2) * v4; }
template<> inline double omega<8,1>(double u, double v) { double v2 = v * v; double v4 = v2 * v2; return 0.013447026066676145 + -1.0516269686752923e-09 * u + (-0.00022648009365806407 + 2.4786189349797326e-11 * u) * v + (4.6005385071969854e-06 + -5.6595863589756684e-13 * u + (-9.9322110199857817e-08 + 1.296787111084816e-14 * u) * v) * v2 + (2.207477230578587e-09 + -2.9875650117833828e-16 * u + (-4.9889965903420901e-11 + 6.9157921515347024e-18 * u) * v + (1.1363652279388493e-12 + -1.6022926068411249e-19 * u + (-2.3949380995435413e-14 + 3.4174991916280297e-21 * u) * v) * v2) * v4; }
template<> inline double omega<8,2>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.020238029360896314 + -3.16900010555155e-09 * u + 1.6503641935804991e-16 * u2 + (-0.00047876629077055698 + 1.0344708713742295e-10 * u + -6.1311683774290491e-18 * u2) * v + (1.2755098155244379e-05 + -3.1085004528166766e-12 * u + 1.9817763011262371e-19 * u2 + (-3.4332369880777869e-07 + 8.9046479350052009e-14 * u + -5.929955027987657e-21 * u2) * v) * v2 + (9.1739225832514903e-09 + -2.470875093941957e-15 * u + 1.6930694246060705e-22 * u2 + (-2.4269904734854989e-10 + 6.7033474197686728e-17 * u + -4.6859311578822256e-24 * u2) * v + (6.3364530639219839e-12 + -1.7816207157382202e-18 * u + 1.2638164838261235e-25 * u2 + (-1.4933476150210949e-13 + 4.2512962699359034e-20 * u + -3.0473670383517657e-27 * u2) * v) * v2) * v4; }
template<> inline double omega<8,3>(double u, double v) { double u2 = u * u; double v2 = v * v; double v4 = v2 * v2; return 0.027074538455918874 + -6.366329240055028e-09 * u + (6.6347138726371738e-16 + -2.5899887028558188e-23 * u) * u2 + (-0.00082615358745242648 + 2.6591508878932403e-10 * u + (-3.145853010338335e-17 + 1.3160739000435209e-24 * u) * u2) * v + (2.7361023869946904e-05 + -9.9394966474793564e-12 * u + (1.2630736601334851e-18 + -5.5334830194186618e-26 * u) * u2 + (-8.875583051213878e-07 + 3.4357838294425397e-13 * u + (-4.5609483796271083e-20 + 2.0618591536880933e-27 * u) * u2) * v) * v2 + (2.7876803473982034e-08 + -1.1218008302703135e-14 * u + (1.5329200946174922e-21 + -7.085733779305315e-29 * u) * u2 + (-8.4974606300087386e-10 + 3.5094351503781681e-16 * u + (-4.8945520968955438e-23 + 2.299864374428438e-30 * u) * u2) * v + (2.5135485567116819e-11 + -1.0573771446003313e-17 * u + (1.4970202530064588e-24 + -7.1222793456298752e-32 * u) * u2 + (-6.5649284381048373e-13 + 2.797071694045094e-19 * u + (-4.0024481079717931e-26 + 1.9214547287321458e-33 * u) * u2) * v) * v2) * v4; }
template<> inline double omega<8,4>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.033956858282643852 + -1.0657992640551349e-08 * u + (1.6670252745297175e-15 + -1.3019553824104611e-22 * u) * u2 + 4.0645825327565944e-30 * u4 + (-0.0012704737749403441 + 5.4298375586824196e-10 * u + (-9.6247791426000922e-17 + 8.0493914875251272e-24 * u) * u2 + -2.6240941524487515e-31 * u4) * v + (5.0409890187223832e-05 + -2.4312302300055643e-11 * u + (4.6239906072118148e-18 + -4.0462347372720098e-25 * u) * u2 + 1.3612421687245953e-32 * u4 + (-1.9192707583887034e-06 + 9.8693431158211243e-13 * u + (-1.960511675526664e-19 + 1.7696774027553613e-26 * u) * u2 + -6.0932233855475119e-34 * u4) * v) * v2 + (6.9512264622097312e-08 + -3.7181117010588571e-14 * u + (7.6041160098051175e-21 + -7.0181796717846223e-28 * u) * u2 + 2.4589772912559127e-35 * u4 + (-2.4070394434672276e-09 + 1.3220446193991829e-15 * u + (-2.7602751638086719e-22 + 2.5899832320376226e-29 * u) * u2 + -9.1973850349728917e-37 * u4) * v + (7.9796560311596872e-11 + -4.4659523612766177e-17 * u + (9.4676354105828963e-24 + -8.9958687064077996e-31 * u) * u2 + 3.2283076617482314e-38 * u4 + (-2.2900009137675096e-12 + 1.2983375871168426e-18 * u + (-2.7822145525216757e-25 + 2.6676646637981374e-32 * u) * u2 + -9.6476290838974729e-40 * u4) * v) * v2) * v4; }
template<> inline double omega<8,5>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.040885295839391086 + -1.6058512798745589e-08 * u + (3.3508256374417514e-15 + -3.9268407515259631e-22 * u) * u2 + (2.452407338013582e-29 + -6.3787271349251216e-37 * u) * u4 + (-0.0018135837717026657 + 9.6632292926978982e-10 * u + (-2.2823864616849956e-16 + 2.8625534089589787e-23 * u) * u2 + (-1.8662042747055426e-30 + 5.0002117515795902e-38 * u) * u4) * v + (8.3959101991121129e-05 + -5.0463006276889892e-11 * u + (1.2777398524239367e-17 + -1.6756648633399439e-24 * u) * u2 + (1.126830577564747e-31 + -3.0887654660146524e-39 * u) * u4 + (-3.6758033943606974e-06 + 2.3561053705861775e-12 * u + (-6.2293914218285041e-19 + 8.4244752487863321e-26 * u) * u2 + (-5.7963737188994945e-33 + 1.6171896595541893e-40 * u) * u4) * v) * v2 + (1.5111191433682134e-07 + -1.0078596565508521e-13 * u + (2.7435219254948867e-20 + -3.7933109403384254e-27 * u) * u2 + (2.6555500693512164e-34 + -7.5126859660501234e-42 * u) * u4 + (-5.8716689704047573e-09 + 4.022637925090707e-15 * u + (-1.1180356647166262e-21 + 1.5716156655614107e-28 * u) * u2 + (-1.1150894761205814e-35 + 3.1898188529147932e-43 * u) * u4) * v + (2.160333038180261e-10 + -1.5085281905170018e-16 * u + (4.2577291172023533e-23 + -6.0610991955731459e-30 * u) * u2 + (4.3459602276482621e-37 + -1.2543125589293998e-44 * u) * u4 + (-6.7587261927338678e-12 + 4.7816057921275127e-18 * u + (-1.3642591458917383e-24 + 1.959808012803829e-31 * u) * u2 + (-1.4161150713137878e-38 + 4.1143122342375957e-46 * u) * u4) * v) * v2) * v4; }
template<> inline double omega<8,6>(double u, double v) { double u2 = u * u; double u4 = u2 * u2; double v2 = v * v; double v4 = v2 * v2; return 0.047860160127022919 + -2.258255765941847e-08 * u + (5.893452665377204e-15 + -9.2118152132823083e-22 * u) * u2 + (8.6314534909413972e-29 + -4.4908326978063562e-36 * u + 1.0010415453327871e-43 * u2) * u4 + (-0.0024573659093425492 + 1.5684946656086881e-09 * u + (-4.6292419937876057e-16 + 7.7406102657791685e-23 * u) * u2 + (-7.5694971244104606e-30 + 4.0563261690809652e-37 * u + -9.2396787364733404e-45 * u2) * u4) * v + (0.00013013328645141266 + -9.3652783333205649e-11 * u + (2.9609844706069168e-17 + -5.1743754588370812e-24 * u) * u2 + (5.2174315577416074e-31 + -2.859584459014862e-38 * u + 6.6270292056421958e-46 * u2) * u4 + (-6.4469256420028762e-06 + 4.9481215703736582e-12 * u + (-1.6331319719219871e-18 + 2.9421867678008347e-25 * u) * u2 + (-3.0345992956755421e-32 + 1.6925197278281919e-39 * u + -3.9771971783987543e-47 * u2) * u4) * v) * v2 + (2.969482826197679e-07 + -2.3719658373229267e-13 * u + (8.0598901844294583e-20 + -1.4843671516899255e-26 * u) * u2 + (1.5575360892836131e-33 + -8.8074710810103531e-41 * u + 2.093095942079841e-48 * u2) * u4 + (-1.2812149298080657e-08 + 1.0514543324522786e-14 * u + (-3.6481455376708318e-21 + 6.830590508825154e-28 * u) * u2 + (-7.2637585206619192e-35 + 4.1530131952920935e-42 * u + -9.9612705695999318e-50 * u2) * u4) * v + (5.1866483166716852e-10 + -4.3392497138742756e-16 * u + (1.5289993663540334e-22 + -2.8992294910981029e-29 * u) * u2 + (3.1156591931186537e-36 + -1.7972219599690599e-43 * u + 4.343514834373485e-51 * u2) * u4 + (-1.7564892580845675e-11 + 1.4889365638943462e-17 * u + (-5.3035924075854948e-24 + 1.014794395797517e-30 * u) * u2 + (-1.0989503454469454e-37 + 6.380970654828389e-45 * u + -1.5509597102505306e-52 * u2) * u4) * v) * v2) * v4; }

}  // namespace
}  // namespace priv

#endif // ECEF2GEO_OMEGA_HPP