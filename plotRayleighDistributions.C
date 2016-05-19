#include "AveragePowerSpectrum.h"

// TString fileNames[] = {"powerSpectraPlots_128_2015-11-23_18-16-33.root", "powerSpectraPlots_130_2015-11-23_18-16-47.root", "powerSpectraPlots_131_2015-11-23_18-16-35.root", "powerSpectraPlots_143_2015-11-23_18-18-07.root", "powerSpectraPlots_142_2015-11-23_18-17-57.root", "powerSpectraPlots_144_2015-11-23_18-18-20.root", "powerSpectraPlots_394_2015-11-24_05-23-02.root", "powerSpectraPlots_132_2015-11-23_18-16-34.root", "powerSpectraPlots_392_2015-11-24_05-20-41.root", "powerSpectraPlots_146_2015-11-23_18-18-27.root", "powerSpectraPlots_395_2015-11-24_05-25-12.root", "powerSpectraPlots_391_2015-11-24_05-20-41.root", "powerSpectraPlots_129_2015-11-23_18-16-57.root", "powerSpectraPlots_133_2015-11-23_18-16-35.root", "powerSpectraPlots_390_2015-11-24_05-20-17.root", "powerSpectraPlots_396_2015-11-24_05-27-43.root", "powerSpectraPlots_134_2015-11-23_18-16-36.root", "powerSpectraPlots_393_2015-11-24_05-31-34.root", "powerSpectraPlots_389_2015-11-24_05-19-12.root", "powerSpectraPlots_294_2015-11-24_00-51-38.root", "powerSpectraPlots_397_2015-11-24_05-27-28.root", "powerSpectraPlots_343_2015-11-24_03-46-58.root", "powerSpectraPlots_386_2015-11-24_05-08-14.root", "powerSpectraPlots_289_2015-11-24_00-08-11.root", "powerSpectraPlots_288_2015-11-24_00-07-43.root", "powerSpectraPlots_293_2015-11-24_00-49-42.root", "powerSpectraPlots_419_2015-11-24_05-52-09.root", "powerSpectraPlots_408_2015-11-24_05-44-24.root", "powerSpectraPlots_291_2015-11-24_00-42-24.root", "powerSpectraPlots_380_2015-11-24_04-58-43.root", "powerSpectraPlots_407_2015-11-24_05-42-17.root", "powerSpectraPlots_418_2015-11-24_05-51-59.root", "powerSpectraPlots_344_2015-11-24_03-47-30.root", "powerSpectraPlots_147_2015-11-23_18-17-33.root", "powerSpectraPlots_342_2015-11-24_03-46-01.root", "powerSpectraPlots_406_2015-11-24_05-40-09.root", "powerSpectraPlots_356_2015-11-24_04-11-59.root", "powerSpectraPlots_352_2015-11-23_16-35-02.root", "powerSpectraPlots_355_2015-11-24_04-11-15.root", "powerSpectraPlots_381_2015-11-24_04-59-11.root", "powerSpectraPlots_352_2015-11-24_04-05-38.root", "powerSpectraPlots_357_2015-11-24_04-13-17.root", "powerSpectraPlots_421_2015-11-24_05-54-28.root", "powerSpectraPlots_417_2015-11-24_05-50-40.root", "powerSpectraPlots_340_2015-11-24_03-39-57.root", "powerSpectraPlots_379_2015-11-24_04-54-54.root", "powerSpectraPlots_345_2015-11-24_03-56-45.root", "powerSpectraPlots_353_2015-11-24_04-06-46.root", "powerSpectraPlots_354_2015-11-24_04-10-13.root", "powerSpectraPlots_319_2015-11-24_02-55-00.root", "powerSpectraPlots_277_2015-11-23_23-18-11.root", "powerSpectraPlots_232_2015-11-23_21-16-10.root", "powerSpectraPlots_292_2015-11-24_00-41-31.root", "powerSpectraPlots_329_2015-11-24_03-13-11.root", "powerSpectraPlots_341_2015-11-24_03-41-17.root", "powerSpectraPlots_196_2015-11-23_19-05-53.root", "powerSpectraPlots_295_2015-11-24_00-55-10.root", "powerSpectraPlots_364_2015-11-24_04-34-51.root", "powerSpectraPlots_382_2015-11-24_05-02-09.root", "powerSpectraPlots_430_2015-11-24_06-05-29.root", "powerSpectraPlots_179_2015-11-23_18-34-25.root", "powerSpectraPlots_318_2015-11-24_02-54-30.root", "powerSpectraPlots_405_2015-11-24_05-40-37.root", "powerSpectraPlots_416_2015-11-24_05-50-38.root", "powerSpectraPlots_282_2015-11-23_23-58-02.root", "powerSpectraPlots_305_2015-11-24_02-12-42.root", "powerSpectraPlots_246_2015-11-23_22-00-02.root", "powerSpectraPlots_431_2015-11-24_06-05-34.root", "powerSpectraPlots_339_2015-11-24_03-38-47.root", "powerSpectraPlots_165_2015-11-23_18-17-12.root", "powerSpectraPlots_178_2015-11-23_18-33-58.root", "powerSpectraPlots_264_2015-11-23_22-32-37.root", "powerSpectraPlots_201_2015-11-23_19-20-45.root", "powerSpectraPlots_409_2015-11-24_05-44-32.root", "powerSpectraPlots_195_2015-11-23_18-59-09.root", "powerSpectraPlots_320_2015-11-24_02-56-46.root", "powerSpectraPlots_245_2015-11-23_21-53-06.root", "powerSpectraPlots_332_2015-11-24_03-20-02.root", "powerSpectraPlots_422_2015-11-24_05-56-38.root", "powerSpectraPlots_333_2015-11-24_03-24-10.root", "powerSpectraPlots_276_2015-11-23_23-06-48.root", "powerSpectraPlots_194_2015-11-23_18-59-20.root", "powerSpectraPlots_231_2015-11-23_21-12-13.root", "powerSpectraPlots_328_2015-11-24_03-10-58.root", "powerSpectraPlots_187_2015-11-23_18-50-54.root", "powerSpectraPlots_346_2015-11-24_03-56-57.root", "powerSpectraPlots_358_2015-11-24_04-14-15.root", "powerSpectraPlots_278_2015-11-23_23-32-13.root", "powerSpectraPlots_188_2015-11-23_18-49-58.root", "powerSpectraPlots_378_2015-11-24_04-55-18.root", "powerSpectraPlots_234_2015-11-23_21-18-01.root", "powerSpectraPlots_275_2015-11-23_22-59-35.root", "powerSpectraPlots_330_2015-11-24_03-13-19.root", "powerSpectraPlots_296_2015-11-24_01-03-08.root", "powerSpectraPlots_365_2015-11-24_04-35-45.root", "powerSpectraPlots_338_2015-11-24_03-38-24.root", "powerSpectraPlots_200_2015-11-23_19-18-18.root", "powerSpectraPlots_331_2015-11-24_03-18-08.root", "powerSpectraPlots_274_2015-11-23_22-56-17.root", "powerSpectraPlots_429_2015-11-24_06-04-39.root", "powerSpectraPlots_366_2015-11-24_04-36-57.root", "powerSpectraPlots_166_2015-11-23_18-17-19.root", "powerSpectraPlots_185_2015-11-23_18-47-16.root", "powerSpectraPlots_309_2015-11-24_02-20-46.root", "powerSpectraPlots_383_2015-11-24_05-04-13.root", "powerSpectraPlots_186_2015-11-23_18-49-20.root", "powerSpectraPlots_398_2015-11-24_05-27-42.root", "powerSpectraPlots_310_2015-11-24_02-22-32.root", "powerSpectraPlots_404_2015-11-24_05-36-55.root", "powerSpectraPlots_281_2015-11-23_23-52-39.root", "powerSpectraPlots_321_2015-11-24_02-57-45.root", "powerSpectraPlots_279_2015-11-23_23-42-47.root", "powerSpectraPlots_248_2015-11-23_22-06-07.root", "powerSpectraPlots_322_2015-11-24_03-00-42.root", "powerSpectraPlots_247_2015-11-23_21-59-10.root", "powerSpectraPlots_164_2015-11-23_18-17-07.root", "powerSpectraPlots_388_2015-11-24_05-18-00.root", "powerSpectraPlots_307_2015-11-24_02-17-44.root", "powerSpectraPlots_308_2015-11-24_02-17-59.root", "powerSpectraPlots_214_2015-11-23_20-15-20.root", "powerSpectraPlots_189_2015-11-23_18-52-14.root", "powerSpectraPlots_362_2015-11-24_04-25-52.root", "powerSpectraPlots_403_2015-11-24_05-36-08.root", "powerSpectraPlots_265_2015-11-23_22-37-54.root", "powerSpectraPlots_387_2015-11-24_05-09-42.root", "powerSpectraPlots_367_2015-11-24_04-37-53.root", "powerSpectraPlots_177_2015-11-23_18-33-54.root", "powerSpectraPlots_297_2015-11-24_01-11-32.root", "powerSpectraPlots_360_2015-11-24_04-16-24.root", "powerSpectraPlots_337_2015-11-24_03-36-54.root", "powerSpectraPlots_174_2015-11-23_18-29-36.root", "powerSpectraPlots_327_2015-11-24_03-11-03.root", "powerSpectraPlots_160_2015-11-23_18-17-08.root", "powerSpectraPlots_317_2015-11-24_02-55-32.root", "powerSpectraPlots_301_2015-11-24_01-54-06.root", "powerSpectraPlots_306_2015-11-24_02-16-34.root", "powerSpectraPlots_176_2015-11-23_18-30-51.root", "powerSpectraPlots_336_2015-11-24_03-35-13.root", "powerSpectraPlots_423_2015-11-24_05-57-15.root", "powerSpectraPlots_361_2015-11-24_04-17-41.root", "powerSpectraPlots_168_2015-11-23_18-17-25.root", "powerSpectraPlots_302_2015-11-24_01-57-20.root", "powerSpectraPlots_350_2015-11-24_04-02-43.root", "powerSpectraPlots_235_2015-11-23_21-20-15.root", "powerSpectraPlots_363_2015-11-24_04-31-53.root", "powerSpectraPlots_304_2015-11-24_02-09-48.root", "powerSpectraPlots_359_2015-11-24_04-16-50.root", "powerSpectraPlots_303_2015-11-24_02-06-20.root", "powerSpectraPlots_167_2015-11-23_18-17-27.root", "powerSpectraPlots_313_2015-11-24_02-33-11.root", "powerSpectraPlots_410_2015-11-24_05-45-38.root", "powerSpectraPlots_300_2015-11-24_01-41-04.root", "powerSpectraPlots_351_2015-11-24_04-04-00.root", "powerSpectraPlots_347_2015-11-24_04-00-51.root", "powerSpectraPlots_298_2015-11-24_01-26-47.root", "powerSpectraPlots_335_2015-11-24_03-29-24.root", "powerSpectraPlots_324_2015-11-24_03-04-42.root", "powerSpectraPlots_163_2015-11-23_18-17-08.root", "powerSpectraPlots_334_2015-11-24_03-27-19.root", "powerSpectraPlots_369_2015-11-24_04-41-47.root", "powerSpectraPlots_280_2015-11-23_23-41-57.root", "powerSpectraPlots_314_2015-11-24_02-37-25.root", "powerSpectraPlots_402_2015-11-24_05-34-52.root", "powerSpectraPlots_249_2015-11-23_22-06-15.root", "powerSpectraPlots_315_2015-11-24_02-49-15.root", "powerSpectraPlots_209_2015-11-23_19-53-05.root", "powerSpectraPlots_311_2015-11-24_02-26-03.root", "powerSpectraPlots_424_2015-11-24_05-58-30.root", "powerSpectraPlots_371_2015-11-24_04-43-31.root", "powerSpectraPlots_427_2015-11-24_06-02-27.root", "powerSpectraPlots_173_2015-11-23_18-24-15.root", "powerSpectraPlots_230_2015-11-23_21-10-09.root", "powerSpectraPlots_299_2015-11-24_01-31-09.root", "powerSpectraPlots_316_2015-11-24_02-54-50.root", "powerSpectraPlots_325_2015-11-24_03-07-06.root", "powerSpectraPlots_348_2015-11-24_04-01-56.root", "powerSpectraPlots_370_2015-11-24_04-42-56.root", "powerSpectraPlots_377_2015-11-24_04-54-11.root", "powerSpectraPlots_428_2015-11-24_06-02-49.root", "powerSpectraPlots_323_2015-11-24_03-02-22.root", "powerSpectraPlots_273_2015-11-23_22-55-14.root", "powerSpectraPlots_215_2015-11-23_20-20-40.root", "powerSpectraPlots_399_2015-11-24_05-28-53.root", "powerSpectraPlots_375_2015-11-24_04-53-50.root", "powerSpectraPlots_349_2015-11-24_04-02-20.root", "powerSpectraPlots_228_2015-11-23_20-54-25.root", "powerSpectraPlots_374_2015-11-24_04-44-47.root", "powerSpectraPlots_266_2015-11-23_22-39-13.root", "powerSpectraPlots_218_2015-11-23_20-30-48.root", "powerSpectraPlots_202_2015-11-23_19-25-26.root", "powerSpectraPlots_272_2015-11-23_22-50-18.root", "powerSpectraPlots_312_2015-11-24_02-33-09.root", "powerSpectraPlots_373_2015-11-24_04-44-03.root", "powerSpectraPlots_400_2015-11-24_05-31-32.root", "powerSpectraPlots_242_2015-11-23_21-36-45.root", "powerSpectraPlots_210_2015-11-23_19-55-47.root", "powerSpectraPlots_184_2015-11-23_18-47-07.root", "powerSpectraPlots_175_2015-11-23_18-29-02.root", "powerSpectraPlots_236_2015-11-23_21-20-56.root", "powerSpectraPlots_415_2015-11-24_05-50-05.root", "powerSpectraPlots_243_2015-11-23_21-50-34.root", "powerSpectraPlots_372_2015-11-24_04-44-51.root", "powerSpectraPlots_211_2015-11-23_20-00-37.root", "powerSpectraPlots_412_2015-11-24_05-47-43.root", "powerSpectraPlots_271_2015-11-23_22-43-52.root", "powerSpectraPlots_385_2015-11-24_05-07-55.root", "powerSpectraPlots_401_2015-11-24_05-33-45.root", "powerSpectraPlots_434_2015-11-24_06-21-59.root", "powerSpectraPlots_238_2015-11-23_21-21-53.root", "powerSpectraPlots_376_2015-11-24_04-54-20.root", "powerSpectraPlots_270_2015-11-23_22-43-33.root", "powerSpectraPlots_326_2015-11-24_03-07-52.root", "powerSpectraPlots_237_2015-11-23_21-22-00.root", "powerSpectraPlots_226_2015-11-23_20-51-36.root", "powerSpectraPlots_207_2015-11-23_19-45-28.root", "powerSpectraPlots_193_2015-11-23_18-54-40.root", "powerSpectraPlots_411_2015-11-24_05-46-32.root", "powerSpectraPlots_213_2015-11-23_20-13-52.root", "powerSpectraPlots_169_2015-11-23_18-17-37.root", "powerSpectraPlots_191_2015-11-23_18-52-43.root", "powerSpectraPlots_425_2015-11-24_05-58-58.root", "powerSpectraPlots_413_2015-11-24_05-49-04.root", "powerSpectraPlots_206_2015-11-23_19-36-20.root", "powerSpectraPlots_414_2015-11-24_05-48-48.root", "powerSpectraPlots_250_2015-11-23_22-08-51.root", "powerSpectraPlots_251_2015-11-23_22-15-39.root", "powerSpectraPlots_192_2015-11-23_18-53-08.root", "powerSpectraPlots_269_2015-11-23_22-40-56.root", "powerSpectraPlots_426_2015-11-24_06-02-22.root", "powerSpectraPlots_241_2015-11-23_21-31-25.root", "powerSpectraPlots_170_2015-11-23_18-17-43.root", "powerSpectraPlots_239_2015-11-23_21-22-55.root", "powerSpectraPlots_254_2015-11-23_22-21-07.root", "powerSpectraPlots_158_2015-11-23_18-17-02.root", "powerSpectraPlots_216_2015-11-23_20-25-44.root", "powerSpectraPlots_224_2015-11-23_20-46-49.root", "powerSpectraPlots_252_2015-11-23_22-18-40.root", "powerSpectraPlots_253_2015-11-23_22-17-53.root", "powerSpectraPlots_432_2015-11-24_06-14-54.root", "powerSpectraPlots_204_2015-11-23_19-30-24.root", "powerSpectraPlots_433_2015-11-24_06-19-22.root", "powerSpectraPlots_197_2015-11-23_19-06-01.root", "powerSpectraPlots_267_2015-11-23_22-39-18.root", "powerSpectraPlots_145_2015-11-23_18-16-51.root", "powerSpectraPlots_152_2015-11-23_17-32-19.root", "powerSpectraPlots_152_2015-11-23_18-00-03.root", "powerSpectraPlots_152_2015-11-23_18-17-07.root", "powerSpectraPlots_255_2015-11-23_22-26-29.root", "powerSpectraPlots_256_2015-11-23_22-28-43.root", "powerSpectraPlots_212_2015-11-23_20-02-06.root", "powerSpectraPlots_171_2015-11-23_18-17-44.root", "powerSpectraPlots_420_2015-11-24_05-54-17.root", "powerSpectraPlots_159_2015-11-23_18-16-58.root", "powerSpectraPlots_290_2015-11-24_00-31-13.root", "powerSpectraPlots_217_2015-11-23_20-29-31.root", "powerSpectraPlots_157_2015-11-23_18-17-01.root", "powerSpectraPlots_208_2015-11-23_19-44-26.root", "powerSpectraPlots_227_2015-11-23_20-53-29.root", "powerSpectraPlots_135_2015-11-23_18-16-46.root", "powerSpectraPlots_141_2015-11-23_18-16-52.root", "powerSpectraPlots_384_2015-11-24_05-05-07.root", "powerSpectraPlots_219_2015-11-23_20-34-09.root", "powerSpectraPlots_149_2015-11-23_18-17-01.root", "powerSpectraPlots_233_2015-11-23_21-16-40.root", "powerSpectraPlots_199_2015-11-23_19-10-47.root", "powerSpectraPlots_284_2015-11-23_23-55-03.root", "powerSpectraPlots_136_2015-11-23_18-16-46.root", "powerSpectraPlots_205_2015-11-23_19-34-40.root", "powerSpectraPlots_222_2015-11-23_20-41-38.root", "powerSpectraPlots_220_2015-11-23_20-37-38.root", "powerSpectraPlots_139_2015-11-23_18-16-52.root", "powerSpectraPlots_140_2015-11-23_18-16-50.root", "powerSpectraPlots_225_2015-11-23_20-45-49.root", "powerSpectraPlots_285_2015-11-24_00-00-50.root", "powerSpectraPlots_190_2015-11-23_18-52-00.root", "powerSpectraPlots_223_2015-11-23_20-44-14.root", "powerSpectraPlots_148_2015-11-23_18-16-57.root", "powerSpectraPlots_268_2015-11-23_22-39-36.root", "powerSpectraPlots_244_2015-11-23_22-03-07.root", "powerSpectraPlots_180_2015-11-23_18-39-11.root", "powerSpectraPlots_161_2015-11-23_18-16-54.root", "powerSpectraPlots_198_2015-11-23_19-08-26.root", "powerSpectraPlots_154_2015-11-23_18-16-53.root", "powerSpectraPlots_151_2015-11-23_18-16-53.root", "powerSpectraPlots_229_2015-11-23_21-08-55.root", "powerSpectraPlots_172_2015-11-23_18-22-25.root", "powerSpectraPlots_221_2015-11-23_20-41-05.root", "powerSpectraPlots_153_2015-11-23_18-16-51.root", "powerSpectraPlots_203_2015-11-23_19-25-15.root", "powerSpectraPlots_138_2015-11-23_18-16-51.root", "powerSpectraPlots_181_2015-11-23_18-39-22.root", "powerSpectraPlots_150_2015-11-23_18-16-49.root", "powerSpectraPlots_162_2015-11-23_18-16-56.root", "powerSpectraPlots_183_2015-11-23_18-43-23.root", "powerSpectraPlots_137_2015-11-23_18-16-51.root", "powerSpectraPlots_155_2015-11-23_18-16-55.root", "powerSpectraPlots_368_2015-11-24_04-42-00.root", "powerSpectraPlots_240_2015-11-23_21-28-23.root", "powerSpectraPlots_283_2015-11-23_23-54-44.root", "powerSpectraPlots_156_2015-11-23_18-16-54.root", "powerSpectraPlots_182_2015-11-23_18-42-03.root", "powerSpectraPlots_286_2015-11-24_00-03-31.root", "powerSpectraPlots_287_2015-11-24_00-05-13.root", "powerSpectraPlots_127_2015-11-23_18-16-33.root", "end"};


const double xSize=750;
const double ySize=625;

void plotRayleighDistributions(){

  // auto f = TFile::Open("powerSpectraPlots_352_2015-11-23_16-35-02.root");
  // TString fileNames[2] = {"powerSpectraPlots_352_2015-11-23_16-35-02.root", "end"};
  TString fileNames[] = {"powerSpectraPlots_270_2015-11-23_22-43-33.root", "end"};
  // TString fileNames[] = {"../powerSpectraPlots_352_2015-11-26_12-43-34.root", "end"};
  // TString fileNames[] = {"../powerSpectraPlots_352_2015-11-26_12-52-16.root", "end"};  

  // TString fileNames[] = {"powerSpectraPlots_161_2015-11-23_18-16-54.root", "powerSpectraPlots_145_2015-11-23_18-16-51.root", "powerSpectraPlots_150_2015-11-23_18-16-49.root", "powerSpectraPlots_135_2015-11-23_18-16-46.root", "powerSpectraPlots_156_2015-11-23_18-16-54.root", "powerSpectraPlots_136_2015-11-23_18-16-46.root", "powerSpectraPlots_162_2015-11-23_18-16-56.root", "powerSpectraPlots_127_2015-11-23_18-16-33.root", "powerSpectraPlots_133_2015-11-23_18-16-35.root", "powerSpectraPlots_155_2015-11-23_18-16-55.root", "powerSpectraPlots_134_2015-11-23_18-16-36.root", "powerSpectraPlots_182_2015-11-23_18-42-03.root", "powerSpectraPlots_180_2015-11-23_18-39-11.root", "powerSpectraPlots_183_2015-11-23_18-43-23.root", "powerSpectraPlots_181_2015-11-23_18-39-22.root", "powerSpectraPlots_139_2015-11-23_18-16-52.root", "powerSpectraPlots_128_2015-11-23_18-16-33.root", "powerSpectraPlots_153_2015-11-23_18-16-51.root", "powerSpectraPlots_138_2015-11-23_18-16-51.root", "powerSpectraPlots_140_2015-11-23_18-16-50.root", "powerSpectraPlots_132_2015-11-23_18-16-34.root", "powerSpectraPlots_137_2015-11-23_18-16-51.root", "powerSpectraPlots_151_2015-11-23_18-16-53.root", "powerSpectraPlots_179_2015-11-23_18-34-25.root", "powerSpectraPlots_141_2015-11-23_18-16-52.root", "powerSpectraPlots_154_2015-11-23_18-16-53.root", "powerSpectraPlots_159_2015-11-23_18-16-58.root", "powerSpectraPlots_198_2015-11-23_19-08-26.root", "powerSpectraPlots_131_2015-11-23_18-16-35.root", "powerSpectraPlots_172_2015-11-23_18-22-25.root", "powerSpectraPlots_148_2015-11-23_18-16-57.root", "powerSpectraPlots_197_2015-11-23_19-06-01.root", "powerSpectraPlots_199_2015-11-23_19-10-47.root", "powerSpectraPlots_190_2015-11-23_18-52-00.root", "powerSpectraPlots_175_2015-11-23_18-29-02.root", "powerSpectraPlots_352_2015-11-23_16-35-02.root", "powerSpectraPlots_152_2015-11-23_18-00-03.root", "end"};


  const Int_t numSamples = 256;
  const Double_t deltaT = 1./2.6;
  Double_t* freqs = FancyFFTs::getFreqArray(numSamples, deltaT);
  const Int_t numFreqs = FancyFFTs::getNumFreqs(numSamples);

  const Double_t maxPow_dB = 50;
  const Double_t minPow_dB = -5;


  const Int_t numRings = 3;
  const Int_t numPhi = 16;  

  TString ringNames[] = {"T", "M", "B"};
  TString polNames[] = {"H", "V"};

  int fileInd = 0;

  while(fileNames[fileInd] != "end"){

    auto f = TFile::Open("powSpecRootFiles/" + fileNames[fileInd]);

    const Int_t ant = 47;
    const Int_t polInd = 1;
    const int run = 270;

    const int ring = ant/numPhi;
    const int phi = ant%numPhi;

    auto c0 = new TCanvas();
    auto grName = TString::Format("grAvePs_%d_%d_%d_%d", polInd, ant, run, run);
    auto gr = (TGraph*) f->Get(grName);
    gr->Draw();
    gr->GetXaxis()->SetNoExponent(1);
    gr->GetYaxis()->SetNoExponent(1);
    auto grTitle = TString::Format("Average Power spectrum - %d%s%s run %d", phi, ringNames[ring].Data(), polNames[polInd].Data(), run);
    grTitle += "; Frequency (MHz); Power Spectra Density (dB)";
    gr->SetTitle(grTitle);
    gr->SetMaximum(maxPow_dB);
    gr->SetMinimum(minPow_dB);

    TH2D* h2 = nullptr;
    Int_t count1 = 0;
    TGraph* grChiSquarePerNDF = new TGraph();
    TGraph* grMeanFromFit = new TGraph();
  
    TCanvas* c1 = new TCanvas();
    TString baseName = "hAvePowSpec";
    for(Int_t freqInd=0; freqInd<numFreqs; freqInd++){
      // if(freqInd!=23) continue;
      // if(freqs[freqInd] > 0 && freqs[freqInd] < 1.3){
      if(freqs[freqInd] > 0.2 && freqs[freqInd] < 0.3){	
	// if(true){      
	TString histName = TString::Format("%s_%d_%d_%d_%d", baseName.Data(), polInd, ant, freqInd, numFreqs);
	auto h = (TH1D*) f->Get(histName);
	h->Rebin(8);
	Int_t nx = h->GetNbinsX();
      
	if(h2==NULL){
	  TString histName2 = TString::Format("%s2_%d_%d_%d_%d", baseName.Data(), polInd, ant, freqInd, numFreqs);
	  auto h2Title = TString::Format("Sqrt(PSD) as a function of frequency - %d%s%s run %d", phi, ringNames[ring].Data(), polNames[polInd].Data(), run);
	  h2Title += ";Frequency (MHz); sqrt(PSD) (some weird units)";
	  h2 = new TH2D(histName2,
			h2Title,
			numFreqs, 1e3*freqs[0], 1e3*numFreqs*(freqs[1]-freqs[0]),
			64, 0, 2048);
	  h2->GetXaxis()->SetNoExponent(1);
	  h2->GetYaxis()->SetNoExponent(1);
	}

	// TCanvas* ct = nullptr;
	TCanvas* ct = new TCanvas();	
	
	TString fitName = "fit_" + histName;
	Double_t lowEdge = h->GetBinLowEdge(1);
	Double_t highEdge = h->GetBinLowEdge(h->GetNbinsX()+1);
	// cout << lowEdge << "\t" << highEdge << endl;

	auto fit = AveragePowerSpectrum::makeRayleighFunction(fitName, lowEdge, highEdge);
	fit->SetParameter(1, 0.75*h->GetMean());
	
	h->GetXaxis()->SetNoExponent(1);
	h->GetYaxis()->SetNoExponent(1);

	h->Fit(fitName, "Q", "0", 0, 1.5*h->GetMean());

	Double_t freq = 1e3*freqs[freqInd];
	for(int binx=1; binx<=nx; binx++){
	  Double_t sqrtPSD = h->GetXaxis()->GetBinCenter(binx);
	  Double_t numEntries =  h->GetBinContent(binx);

	  Double_t fitVal = fit->Eval(sqrtPSD);
	  h2->Fill(freq, sqrtPSD, numEntries);
	}

	

	// TCanvas* ct = nullptr;

	  
	// auto h = (TH1D*) f->Get(TString::Format("hAvePowSpecs_0_23_%d_%d", freqInd, numFreqs));
	TString opt = "e";
	Double_t chiSq = fit->GetChisquare();
	// std::cout << freqInd << "\t" << chiSq << "\t" << fit->GetNDF() << std::endl;
	chiSq /= fit->GetNDF();
	grChiSquarePerNDF->SetPoint(grChiSquarePerNDF->GetN(), freqs[freqInd]*1e3, chiSq);


	// TGraph* grTest = new TGraph();
	// Double_t mean = h->GetMean();
	// Int_t counter=0;
	// for(int binx=nx; binx>=2; binx--){
	//   TString name99 = TString::Format("clone%d", binx);	  
	//   auto hClone = (TH1D*) h->Clone(name99);
	//   auto f2 = (TF1*) fit->Clone();
	//   for(int binxx=nx; binxx>=binx; binxx--){
	//     hClone->SetBinContent(binxx,0);
	//     hClone->SetBinError(binxx,0);
	//   }
	//   f2->SetParameter(1, mean);
	//   TString name = TString::Format("fit%d", binx);
	//   f2->SetName(name);
	//   hClone->Fit(name, "Q", "0");//, 0, 1.5*h->GetMean());
	//   if(f2->GetNDF() > 0){
	//     grTest->SetPoint(grTest->GetN(), counter, f2->GetChisquare() / f2->GetNDF());
	//     counter++;
	//   }
	//   new TCanvas();
	//   hClone->Draw();
	//   f2->Draw("same");
	// }
	// new TCanvas();
	// grTest->Draw();
	
	
	if(ct!=NULL){
	  ct->Clear();
	  h->Draw(opt);
	  fit->Draw("lsame");
	  //	  gPad->SetLogy(1);
	}
      }
    }
    delete c1;
    auto c2 = new TCanvas();
    h2->Draw("colz");
    gPad->SetLogz(1);
              
    c0->cd();
    grChiSquarePerNDF->SetLineColor(kRed);
    grChiSquarePerNDF->Draw("lsame");
    
    TLegend* l = new TLegend(0.8, 0.8, 1, 1);
    l->AddEntry(gr, "PSD", "l");
    l->AddEntry(grChiSquarePerNDF, "#chi^{2}/NDF Rayleigh Fit", "l");
    l->Draw();
    // f->Close();    
    fileInd++;
  }

  
}
