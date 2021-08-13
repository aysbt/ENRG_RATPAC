// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__RATDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "RAT/DS/Root.hh"
#include "RAT/DS/EV.hh"
#include "RAT/DS/MC.hh"
#include "RAT/DS/MCParticle.hh"
#include "RAT/DS/MCPhoton.hh"
#include "RAT/DS/MCPMT.hh"
#include "RAT/DS/MCSummary.hh"
#include "RAT/DS/PMT.hh"
#include "RAT/DS/RunStore.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/PosFit.hh"
#include "RAT/DS/PMTInfo.hh"
#include "RAT/DS/MCTrack.hh"
#include "RAT/DS/MCTrackStep.hh"
#include "RAT/DS/Calib.hh"
#include "RAT/DS/Centroid.hh"
#include "RAT/DS/PathFit.hh"
#include "RAT/DS/BonsaiFit.hh"
#include "RAT/DSReader.hh"
#include "RAT/DSWriter.hh"
#include "RAT/TrackNav.hh"
#include "RAT/TrackNode.hh"
#include "RAT/TrackCursor.hh"
#include "RAT/DB.hh"
#include "RAT/DBLink.hh"
#include "RAT/DBTextLoader.hh"
#include "RAT/DBJsonLoader.hh"
#include "RAT/HTTPDownloader.hh"
#include "RAT/Log.hh"
#include "RAT/ObjInt.hh"
#include "RAT/ObjDbl.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *pairlEstringcOstringgR_Dictionary();
   static void pairlEstringcOstringgR_TClassManip(TClass*);
   static void *new_pairlEstringcOstringgR(void *p = 0);
   static void *newArray_pairlEstringcOstringgR(Long_t size, void *p);
   static void delete_pairlEstringcOstringgR(void *p);
   static void deleteArray_pairlEstringcOstringgR(void *p);
   static void destruct_pairlEstringcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const pair<string,string>*)
   {
      pair<string,string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(pair<string,string>));
      static ::ROOT::TGenericClassInfo 
         instance("pair<string,string>", "string", 208,
                  typeid(pair<string,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairlEstringcOstringgR_Dictionary, isa_proxy, 4,
                  sizeof(pair<string,string>) );
      instance.SetNew(&new_pairlEstringcOstringgR);
      instance.SetNewArray(&newArray_pairlEstringcOstringgR);
      instance.SetDelete(&delete_pairlEstringcOstringgR);
      instance.SetDeleteArray(&deleteArray_pairlEstringcOstringgR);
      instance.SetDestructor(&destruct_pairlEstringcOstringgR);
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const pair<string,string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairlEstringcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const pair<string,string>*)0x0)->GetClass();
      pairlEstringcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void pairlEstringcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLMCParticle(void *p = 0);
   static void *newArray_RATcLcLDScLcLMCParticle(Long_t size, void *p);
   static void delete_RATcLcLDScLcLMCParticle(void *p);
   static void deleteArray_RATcLcLDScLcLMCParticle(void *p);
   static void destruct_RATcLcLDScLcLMCParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::MCParticle*)
   {
      ::RAT::DS::MCParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::MCParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::MCParticle", ::RAT::DS::MCParticle::Class_Version(), "include/RAT/DS/MCParticle.hh", 21,
                  typeid(::RAT::DS::MCParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::MCParticle::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::MCParticle) );
      instance.SetNew(&new_RATcLcLDScLcLMCParticle);
      instance.SetNewArray(&newArray_RATcLcLDScLcLMCParticle);
      instance.SetDelete(&delete_RATcLcLDScLcLMCParticle);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLMCParticle);
      instance.SetDestructor(&destruct_RATcLcLDScLcLMCParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::MCParticle*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::MCParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::MCParticle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLMCTrackStep(void *p = 0);
   static void *newArray_RATcLcLDScLcLMCTrackStep(Long_t size, void *p);
   static void delete_RATcLcLDScLcLMCTrackStep(void *p);
   static void deleteArray_RATcLcLDScLcLMCTrackStep(void *p);
   static void destruct_RATcLcLDScLcLMCTrackStep(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::MCTrackStep*)
   {
      ::RAT::DS::MCTrackStep *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::MCTrackStep >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::MCTrackStep", ::RAT::DS::MCTrackStep::Class_Version(), "include/RAT/DS/MCTrackStep.hh", 28,
                  typeid(::RAT::DS::MCTrackStep), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::MCTrackStep::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::MCTrackStep) );
      instance.SetNew(&new_RATcLcLDScLcLMCTrackStep);
      instance.SetNewArray(&newArray_RATcLcLDScLcLMCTrackStep);
      instance.SetDelete(&delete_RATcLcLDScLcLMCTrackStep);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLMCTrackStep);
      instance.SetDestructor(&destruct_RATcLcLDScLcLMCTrackStep);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::MCTrackStep*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::MCTrackStep*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::MCTrackStep*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLMCTrack(void *p = 0);
   static void *newArray_RATcLcLDScLcLMCTrack(Long_t size, void *p);
   static void delete_RATcLcLDScLcLMCTrack(void *p);
   static void deleteArray_RATcLcLDScLcLMCTrack(void *p);
   static void destruct_RATcLcLDScLcLMCTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::MCTrack*)
   {
      ::RAT::DS::MCTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::MCTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::MCTrack", ::RAT::DS::MCTrack::Class_Version(), "include/RAT/DS/MCTrack.hh", 21,
                  typeid(::RAT::DS::MCTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::MCTrack::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::MCTrack) );
      instance.SetNew(&new_RATcLcLDScLcLMCTrack);
      instance.SetNewArray(&newArray_RATcLcLDScLcLMCTrack);
      instance.SetDelete(&delete_RATcLcLDScLcLMCTrack);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLMCTrack);
      instance.SetDestructor(&destruct_RATcLcLDScLcLMCTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::MCTrack*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::MCTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::MCTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLMCSummary(void *p = 0);
   static void *newArray_RATcLcLDScLcLMCSummary(Long_t size, void *p);
   static void delete_RATcLcLDScLcLMCSummary(void *p);
   static void deleteArray_RATcLcLDScLcLMCSummary(void *p);
   static void destruct_RATcLcLDScLcLMCSummary(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::MCSummary*)
   {
      ::RAT::DS::MCSummary *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::MCSummary >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::MCSummary", ::RAT::DS::MCSummary::Class_Version(), "include/RAT/DS/MCSummary.hh", 24,
                  typeid(::RAT::DS::MCSummary), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::MCSummary::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::MCSummary) );
      instance.SetNew(&new_RATcLcLDScLcLMCSummary);
      instance.SetNewArray(&newArray_RATcLcLDScLcLMCSummary);
      instance.SetDelete(&delete_RATcLcLDScLcLMCSummary);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLMCSummary);
      instance.SetDestructor(&destruct_RATcLcLDScLcLMCSummary);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::MCSummary*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::MCSummary*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::MCSummary*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLMCPhoton(void *p = 0);
   static void *newArray_RATcLcLDScLcLMCPhoton(Long_t size, void *p);
   static void delete_RATcLcLDScLcLMCPhoton(void *p);
   static void deleteArray_RATcLcLDScLcLMCPhoton(void *p);
   static void destruct_RATcLcLDScLcLMCPhoton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::MCPhoton*)
   {
      ::RAT::DS::MCPhoton *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::MCPhoton >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::MCPhoton", ::RAT::DS::MCPhoton::Class_Version(), "include/RAT/DS/MCPhoton.hh", 25,
                  typeid(::RAT::DS::MCPhoton), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::MCPhoton::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::MCPhoton) );
      instance.SetNew(&new_RATcLcLDScLcLMCPhoton);
      instance.SetNewArray(&newArray_RATcLcLDScLcLMCPhoton);
      instance.SetDelete(&delete_RATcLcLDScLcLMCPhoton);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLMCPhoton);
      instance.SetDestructor(&destruct_RATcLcLDScLcLMCPhoton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::MCPhoton*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::MCPhoton*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::MCPhoton*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *jsoncLcLValue_Dictionary();
   static void jsoncLcLValue_TClassManip(TClass*);
   static void *new_jsoncLcLValue(void *p = 0);
   static void *newArray_jsoncLcLValue(Long_t size, void *p);
   static void delete_jsoncLcLValue(void *p);
   static void deleteArray_jsoncLcLValue(void *p);
   static void destruct_jsoncLcLValue(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::json::Value*)
   {
      ::json::Value *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::json::Value));
      static ::ROOT::TGenericClassInfo 
         instance("json::Value", "RAT/json.hh", 67,
                  typeid(::json::Value), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &jsoncLcLValue_Dictionary, isa_proxy, 4,
                  sizeof(::json::Value) );
      instance.SetNew(&new_jsoncLcLValue);
      instance.SetNewArray(&newArray_jsoncLcLValue);
      instance.SetDelete(&delete_jsoncLcLValue);
      instance.SetDeleteArray(&deleteArray_jsoncLcLValue);
      instance.SetDestructor(&destruct_jsoncLcLValue);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::json::Value*)
   {
      return GenerateInitInstanceLocal((::json::Value*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::json::Value*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *jsoncLcLValue_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::json::Value*)0x0)->GetClass();
      jsoncLcLValue_TClassManip(theClass);
   return theClass;
   }

   static void jsoncLcLValue_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLLog_Dictionary();
   static void RATcLcLLog_TClassManip(TClass*);
   static void delete_RATcLcLLog(void *p);
   static void deleteArray_RATcLcLLog(void *p);
   static void destruct_RATcLcLLog(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::Log*)
   {
      ::RAT::Log *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::Log));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::Log", "RAT/Log.hh", 91,
                  typeid(::RAT::Log), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLLog_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::Log) );
      instance.SetDelete(&delete_RATcLcLLog);
      instance.SetDeleteArray(&deleteArray_RATcLcLLog);
      instance.SetDestructor(&destruct_RATcLcLLog);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::Log*)
   {
      return GenerateInitInstanceLocal((::RAT::Log*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::Log*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLLog_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::Log*)0x0)->GetClass();
      RATcLcLLog_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLLog_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLMCPMT(void *p = 0);
   static void *newArray_RATcLcLDScLcLMCPMT(Long_t size, void *p);
   static void delete_RATcLcLDScLcLMCPMT(void *p);
   static void deleteArray_RATcLcLDScLcLMCPMT(void *p);
   static void destruct_RATcLcLDScLcLMCPMT(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::MCPMT*)
   {
      ::RAT::DS::MCPMT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::MCPMT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::MCPMT", ::RAT::DS::MCPMT::Class_Version(), "include/RAT/DS/MCPMT.hh", 21,
                  typeid(::RAT::DS::MCPMT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::MCPMT::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::MCPMT) );
      instance.SetNew(&new_RATcLcLDScLcLMCPMT);
      instance.SetNewArray(&newArray_RATcLcLDScLcLMCPMT);
      instance.SetDelete(&delete_RATcLcLDScLcLMCPMT);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLMCPMT);
      instance.SetDestructor(&destruct_RATcLcLDScLcLMCPMT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::MCPMT*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::MCPMT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::MCPMT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLMC(void *p = 0);
   static void *newArray_RATcLcLDScLcLMC(Long_t size, void *p);
   static void delete_RATcLcLDScLcLMC(void *p);
   static void deleteArray_RATcLcLDScLcLMC(void *p);
   static void destruct_RATcLcLDScLcLMC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::MC*)
   {
      ::RAT::DS::MC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::MC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::MC", ::RAT::DS::MC::Class_Version(), "include/RAT/DS/MC.hh", 28,
                  typeid(::RAT::DS::MC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::MC::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::MC) );
      instance.SetNew(&new_RATcLcLDScLcLMC);
      instance.SetNewArray(&newArray_RATcLcLDScLcLMC);
      instance.SetDelete(&delete_RATcLcLDScLcLMC);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLMC);
      instance.SetDestructor(&destruct_RATcLcLDScLcLMC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::MC*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::MC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::MC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLPMT(void *p = 0);
   static void *newArray_RATcLcLDScLcLPMT(Long_t size, void *p);
   static void delete_RATcLcLDScLcLPMT(void *p);
   static void deleteArray_RATcLcLDScLcLPMT(void *p);
   static void destruct_RATcLcLDScLcLPMT(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::PMT*)
   {
      ::RAT::DS::PMT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::PMT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::PMT", ::RAT::DS::PMT::Class_Version(), "include/RAT/DS/PMT.hh", 16,
                  typeid(::RAT::DS::PMT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::PMT::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::PMT) );
      instance.SetNew(&new_RATcLcLDScLcLPMT);
      instance.SetNewArray(&newArray_RATcLcLDScLcLPMT);
      instance.SetDelete(&delete_RATcLcLDScLcLPMT);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLPMT);
      instance.SetDestructor(&destruct_RATcLcLDScLcLPMT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::PMT*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::PMT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::PMT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RATcLcLDScLcLPosFit(void *p);
   static void deleteArray_RATcLcLDScLcLPosFit(void *p);
   static void destruct_RATcLcLDScLcLPosFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::PosFit*)
   {
      ::RAT::DS::PosFit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::PosFit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::PosFit", ::RAT::DS::PosFit::Class_Version(), "include/RAT/DS/PosFit.hh", 21,
                  typeid(::RAT::DS::PosFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::PosFit::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::PosFit) );
      instance.SetDelete(&delete_RATcLcLDScLcLPosFit);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLPosFit);
      instance.SetDestructor(&destruct_RATcLcLDScLcLPosFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::PosFit*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::PosFit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::PosFit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLCentroid(void *p = 0);
   static void *newArray_RATcLcLDScLcLCentroid(Long_t size, void *p);
   static void delete_RATcLcLDScLcLCentroid(void *p);
   static void deleteArray_RATcLcLDScLcLCentroid(void *p);
   static void destruct_RATcLcLDScLcLCentroid(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::Centroid*)
   {
      ::RAT::DS::Centroid *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::Centroid >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::Centroid", ::RAT::DS::Centroid::Class_Version(), "include/RAT/DS/Centroid.hh", 14,
                  typeid(::RAT::DS::Centroid), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::Centroid::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::Centroid) );
      instance.SetNew(&new_RATcLcLDScLcLCentroid);
      instance.SetNewArray(&newArray_RATcLcLDScLcLCentroid);
      instance.SetDelete(&delete_RATcLcLDScLcLCentroid);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLCentroid);
      instance.SetDestructor(&destruct_RATcLcLDScLcLCentroid);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::Centroid*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::Centroid*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::Centroid*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLBonsaiFit(void *p = 0);
   static void *newArray_RATcLcLDScLcLBonsaiFit(Long_t size, void *p);
   static void delete_RATcLcLDScLcLBonsaiFit(void *p);
   static void deleteArray_RATcLcLDScLcLBonsaiFit(void *p);
   static void destruct_RATcLcLDScLcLBonsaiFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::BonsaiFit*)
   {
      ::RAT::DS::BonsaiFit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::BonsaiFit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::BonsaiFit", ::RAT::DS::BonsaiFit::Class_Version(), "include/RAT/DS/BonsaiFit.hh", 17,
                  typeid(::RAT::DS::BonsaiFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::BonsaiFit::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::BonsaiFit) );
      instance.SetNew(&new_RATcLcLDScLcLBonsaiFit);
      instance.SetNewArray(&newArray_RATcLcLDScLcLBonsaiFit);
      instance.SetDelete(&delete_RATcLcLDScLcLBonsaiFit);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLBonsaiFit);
      instance.SetDestructor(&destruct_RATcLcLDScLcLBonsaiFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::BonsaiFit*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::BonsaiFit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::BonsaiFit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLPathFit(void *p = 0);
   static void *newArray_RATcLcLDScLcLPathFit(Long_t size, void *p);
   static void delete_RATcLcLDScLcLPathFit(void *p);
   static void deleteArray_RATcLcLDScLcLPathFit(void *p);
   static void destruct_RATcLcLDScLcLPathFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::PathFit*)
   {
      ::RAT::DS::PathFit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::PathFit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::PathFit", ::RAT::DS::PathFit::Class_Version(), "include/RAT/DS/PathFit.hh", 14,
                  typeid(::RAT::DS::PathFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::PathFit::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::PathFit) );
      instance.SetNew(&new_RATcLcLDScLcLPathFit);
      instance.SetNewArray(&newArray_RATcLcLDScLcLPathFit);
      instance.SetDelete(&delete_RATcLcLDScLcLPathFit);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLPathFit);
      instance.SetDestructor(&destruct_RATcLcLDScLcLPathFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::PathFit*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::PathFit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::PathFit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLEV(void *p = 0);
   static void *newArray_RATcLcLDScLcLEV(Long_t size, void *p);
   static void delete_RATcLcLDScLcLEV(void *p);
   static void deleteArray_RATcLcLDScLcLEV(void *p);
   static void destruct_RATcLcLDScLcLEV(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::EV*)
   {
      ::RAT::DS::EV *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::EV >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::EV", ::RAT::DS::EV::Class_Version(), "include/RAT/DS/EV.hh", 30,
                  typeid(::RAT::DS::EV), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::EV::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::EV) );
      instance.SetNew(&new_RATcLcLDScLcLEV);
      instance.SetNewArray(&newArray_RATcLcLDScLcLEV);
      instance.SetDelete(&delete_RATcLcLDScLcLEV);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLEV);
      instance.SetDestructor(&destruct_RATcLcLDScLcLEV);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::EV*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::EV*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::EV*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLCalib(void *p = 0);
   static void *newArray_RATcLcLDScLcLCalib(Long_t size, void *p);
   static void delete_RATcLcLDScLcLCalib(void *p);
   static void deleteArray_RATcLcLDScLcLCalib(void *p);
   static void destruct_RATcLcLDScLcLCalib(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::Calib*)
   {
      ::RAT::DS::Calib *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::Calib >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::Calib", ::RAT::DS::Calib::Class_Version(), "include/RAT/DS/Calib.hh", 29,
                  typeid(::RAT::DS::Calib), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::Calib::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::Calib) );
      instance.SetNew(&new_RATcLcLDScLcLCalib);
      instance.SetNewArray(&newArray_RATcLcLDScLcLCalib);
      instance.SetDelete(&delete_RATcLcLDScLcLCalib);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLCalib);
      instance.SetDestructor(&destruct_RATcLcLDScLcLCalib);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::Calib*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::Calib*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::Calib*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLRoot(void *p = 0);
   static void *newArray_RATcLcLDScLcLRoot(Long_t size, void *p);
   static void delete_RATcLcLDScLcLRoot(void *p);
   static void deleteArray_RATcLcLDScLcLRoot(void *p);
   static void destruct_RATcLcLDScLcLRoot(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::Root*)
   {
      ::RAT::DS::Root *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::Root >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::Root", ::RAT::DS::Root::Class_Version(), "include/RAT/DS/Root.hh", 38,
                  typeid(::RAT::DS::Root), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::Root::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::Root) );
      instance.SetNew(&new_RATcLcLDScLcLRoot);
      instance.SetNewArray(&newArray_RATcLcLDScLcLRoot);
      instance.SetDelete(&delete_RATcLcLDScLcLRoot);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLRoot);
      instance.SetDestructor(&destruct_RATcLcLDScLcLRoot);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::Root*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::Root*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::Root*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLPMTInfo(void *p = 0);
   static void *newArray_RATcLcLDScLcLPMTInfo(Long_t size, void *p);
   static void delete_RATcLcLDScLcLPMTInfo(void *p);
   static void deleteArray_RATcLcLDScLcLPMTInfo(void *p);
   static void destruct_RATcLcLDScLcLPMTInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::PMTInfo*)
   {
      ::RAT::DS::PMTInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::PMTInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::PMTInfo", ::RAT::DS::PMTInfo::Class_Version(), "include/RAT/DS/PMTInfo.hh", 18,
                  typeid(::RAT::DS::PMTInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::PMTInfo::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::PMTInfo) );
      instance.SetNew(&new_RATcLcLDScLcLPMTInfo);
      instance.SetNewArray(&newArray_RATcLcLDScLcLPMTInfo);
      instance.SetDelete(&delete_RATcLcLDScLcLPMTInfo);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLPMTInfo);
      instance.SetDestructor(&destruct_RATcLcLDScLcLPMTInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::PMTInfo*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::PMTInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::PMTInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLRun(void *p = 0);
   static void *newArray_RATcLcLDScLcLRun(Long_t size, void *p);
   static void delete_RATcLcLDScLcLRun(void *p);
   static void deleteArray_RATcLcLDScLcLRun(void *p);
   static void destruct_RATcLcLDScLcLRun(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::Run*)
   {
      ::RAT::DS::Run *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::Run >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::Run", ::RAT::DS::Run::Class_Version(), "include/RAT/DS/Run.hh", 19,
                  typeid(::RAT::DS::Run), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::Run::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::Run) );
      instance.SetNew(&new_RATcLcLDScLcLRun);
      instance.SetNewArray(&newArray_RATcLcLDScLcLRun);
      instance.SetDelete(&delete_RATcLcLDScLcLRun);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLRun);
      instance.SetDestructor(&destruct_RATcLcLDScLcLRun);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::Run*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::Run*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::Run*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLDScLcLRunStore(void *p = 0);
   static void *newArray_RATcLcLDScLcLRunStore(Long_t size, void *p);
   static void delete_RATcLcLDScLcLRunStore(void *p);
   static void deleteArray_RATcLcLDScLcLRunStore(void *p);
   static void destruct_RATcLcLDScLcLRunStore(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DS::RunStore*)
   {
      ::RAT::DS::RunStore *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DS::RunStore >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DS::RunStore", ::RAT::DS::RunStore::Class_Version(), "include/RAT/DS/RunStore.hh", 78,
                  typeid(::RAT::DS::RunStore), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DS::RunStore::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DS::RunStore) );
      instance.SetNew(&new_RATcLcLDScLcLRunStore);
      instance.SetNewArray(&newArray_RATcLcLDScLcLRunStore);
      instance.SetDelete(&delete_RATcLcLDScLcLRunStore);
      instance.SetDeleteArray(&deleteArray_RATcLcLDScLcLRunStore);
      instance.SetDestructor(&destruct_RATcLcLDScLcLRunStore);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DS::RunStore*)
   {
      return GenerateInitInstanceLocal((::RAT::DS::RunStore*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DS::RunStore*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RATcLcLDSReader(void *p);
   static void deleteArray_RATcLcLDSReader(void *p);
   static void destruct_RATcLcLDSReader(void *p);
   static void streamer_RATcLcLDSReader(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DSReader*)
   {
      ::RAT::DSReader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DSReader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DSReader", ::RAT::DSReader::Class_Version(), "RAT/DSReader.hh", 13,
                  typeid(::RAT::DSReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DSReader::Dictionary, isa_proxy, 16,
                  sizeof(::RAT::DSReader) );
      instance.SetDelete(&delete_RATcLcLDSReader);
      instance.SetDeleteArray(&deleteArray_RATcLcLDSReader);
      instance.SetDestructor(&destruct_RATcLcLDSReader);
      instance.SetStreamerFunc(&streamer_RATcLcLDSReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DSReader*)
   {
      return GenerateInitInstanceLocal((::RAT::DSReader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DSReader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RATcLcLDSWriter(void *p);
   static void deleteArray_RATcLcLDSWriter(void *p);
   static void destruct_RATcLcLDSWriter(void *p);
   static void streamer_RATcLcLDSWriter(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DSWriter*)
   {
      ::RAT::DSWriter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::DSWriter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DSWriter", ::RAT::DSWriter::Class_Version(), "RAT/DSWriter.hh", 14,
                  typeid(::RAT::DSWriter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::DSWriter::Dictionary, isa_proxy, 16,
                  sizeof(::RAT::DSWriter) );
      instance.SetDelete(&delete_RATcLcLDSWriter);
      instance.SetDeleteArray(&deleteArray_RATcLcLDSWriter);
      instance.SetDestructor(&destruct_RATcLcLDSWriter);
      instance.SetStreamerFunc(&streamer_RATcLcLDSWriter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DSWriter*)
   {
      return GenerateInitInstanceLocal((::RAT::DSWriter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DSWriter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLTrackNode(void *p = 0);
   static void *newArray_RATcLcLTrackNode(Long_t size, void *p);
   static void delete_RATcLcLTrackNode(void *p);
   static void deleteArray_RATcLcLTrackNode(void *p);
   static void destruct_RATcLcLTrackNode(void *p);
   static void streamer_RATcLcLTrackNode(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::TrackNode*)
   {
      ::RAT::TrackNode *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::TrackNode >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::TrackNode", ::RAT::TrackNode::Class_Version(), "RAT/TrackNode.hh", 10,
                  typeid(::RAT::TrackNode), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::TrackNode::Dictionary, isa_proxy, 16,
                  sizeof(::RAT::TrackNode) );
      instance.SetNew(&new_RATcLcLTrackNode);
      instance.SetNewArray(&newArray_RATcLcLTrackNode);
      instance.SetDelete(&delete_RATcLcLTrackNode);
      instance.SetDeleteArray(&deleteArray_RATcLcLTrackNode);
      instance.SetDestructor(&destruct_RATcLcLTrackNode);
      instance.SetStreamerFunc(&streamer_RATcLcLTrackNode);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::TrackNode*)
   {
      return GenerateInitInstanceLocal((::RAT::TrackNode*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::TrackNode*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLTrackCursor_Dictionary();
   static void RATcLcLTrackCursor_TClassManip(TClass*);
   static void delete_RATcLcLTrackCursor(void *p);
   static void deleteArray_RATcLcLTrackCursor(void *p);
   static void destruct_RATcLcLTrackCursor(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::TrackCursor*)
   {
      ::RAT::TrackCursor *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::TrackCursor));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::TrackCursor", "RAT/TrackCursor.hh", 16,
                  typeid(::RAT::TrackCursor), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLTrackCursor_Dictionary, isa_proxy, 0,
                  sizeof(::RAT::TrackCursor) );
      instance.SetDelete(&delete_RATcLcLTrackCursor);
      instance.SetDeleteArray(&deleteArray_RATcLcLTrackCursor);
      instance.SetDestructor(&destruct_RATcLcLTrackCursor);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::TrackCursor*)
   {
      return GenerateInitInstanceLocal((::RAT::TrackCursor*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::TrackCursor*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLTrackCursor_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::TrackCursor*)0x0)->GetClass();
      RATcLcLTrackCursor_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLTrackCursor_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLTrackNav_Dictionary();
   static void RATcLcLTrackNav_TClassManip(TClass*);
   static void delete_RATcLcLTrackNav(void *p);
   static void deleteArray_RATcLcLTrackNav(void *p);
   static void destruct_RATcLcLTrackNav(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::TrackNav*)
   {
      ::RAT::TrackNav *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::TrackNav));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::TrackNav", "RAT/TrackNav.hh", 11,
                  typeid(::RAT::TrackNav), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLTrackNav_Dictionary, isa_proxy, 0,
                  sizeof(::RAT::TrackNav) );
      instance.SetDelete(&delete_RATcLcLTrackNav);
      instance.SetDeleteArray(&deleteArray_RATcLcLTrackNav);
      instance.SetDestructor(&destruct_RATcLcLTrackNav);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::TrackNav*)
   {
      return GenerateInitInstanceLocal((::RAT::TrackNav*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::TrackNav*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLTrackNav_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::TrackNav*)0x0)->GetClass();
      RATcLcLTrackNav_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLTrackNav_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLHTTPDownloader_Dictionary();
   static void RATcLcLHTTPDownloader_TClassManip(TClass*);
   static void *new_RATcLcLHTTPDownloader(void *p = 0);
   static void *newArray_RATcLcLHTTPDownloader(Long_t size, void *p);
   static void delete_RATcLcLHTTPDownloader(void *p);
   static void deleteArray_RATcLcLHTTPDownloader(void *p);
   static void destruct_RATcLcLHTTPDownloader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::HTTPDownloader*)
   {
      ::RAT::HTTPDownloader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::HTTPDownloader));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::HTTPDownloader", "RAT/HTTPDownloader.hh", 20,
                  typeid(::RAT::HTTPDownloader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLHTTPDownloader_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::HTTPDownloader) );
      instance.SetNew(&new_RATcLcLHTTPDownloader);
      instance.SetNewArray(&newArray_RATcLcLHTTPDownloader);
      instance.SetDelete(&delete_RATcLcLHTTPDownloader);
      instance.SetDeleteArray(&deleteArray_RATcLcLHTTPDownloader);
      instance.SetDestructor(&destruct_RATcLcLHTTPDownloader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::HTTPDownloader*)
   {
      return GenerateInitInstanceLocal((::RAT::HTTPDownloader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::HTTPDownloader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLHTTPDownloader_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::HTTPDownloader*)0x0)->GetClass();
      RATcLcLHTTPDownloader_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLHTTPDownloader_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLDBFieldCallback_Dictionary();
   static void RATcLcLDBFieldCallback_TClassManip(TClass*);
   static void delete_RATcLcLDBFieldCallback(void *p);
   static void deleteArray_RATcLcLDBFieldCallback(void *p);
   static void destruct_RATcLcLDBFieldCallback(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DBFieldCallback*)
   {
      ::RAT::DBFieldCallback *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::DBFieldCallback));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DBFieldCallback", "RAT/DBFieldCallback.hh", 8,
                  typeid(::RAT::DBFieldCallback), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLDBFieldCallback_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DBFieldCallback) );
      instance.SetDelete(&delete_RATcLcLDBFieldCallback);
      instance.SetDeleteArray(&deleteArray_RATcLcLDBFieldCallback);
      instance.SetDestructor(&destruct_RATcLcLDBFieldCallback);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DBFieldCallback*)
   {
      return GenerateInitInstanceLocal((::RAT::DBFieldCallback*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DBFieldCallback*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLDBFieldCallback_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::DBFieldCallback*)0x0)->GetClass();
      RATcLcLDBFieldCallback_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLDBFieldCallback_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLDBTable_Dictionary();
   static void RATcLcLDBTable_TClassManip(TClass*);
   static void *new_RATcLcLDBTable(void *p = 0);
   static void *newArray_RATcLcLDBTable(Long_t size, void *p);
   static void delete_RATcLcLDBTable(void *p);
   static void deleteArray_RATcLcLDBTable(void *p);
   static void destruct_RATcLcLDBTable(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DBTable*)
   {
      ::RAT::DBTable *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::DBTable));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DBTable", "RAT/DBTable.hh", 22,
                  typeid(::RAT::DBTable), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLDBTable_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DBTable) );
      instance.SetNew(&new_RATcLcLDBTable);
      instance.SetNewArray(&newArray_RATcLcLDBTable);
      instance.SetDelete(&delete_RATcLcLDBTable);
      instance.SetDeleteArray(&deleteArray_RATcLcLDBTable);
      instance.SetDestructor(&destruct_RATcLcLDBTable);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DBTable*)
   {
      return GenerateInitInstanceLocal((::RAT::DBTable*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DBTable*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLDBTable_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::DBTable*)0x0)->GetClass();
      RATcLcLDBTable_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLDBTable_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLDBLink_Dictionary();
   static void RATcLcLDBLink_TClassManip(TClass*);
   static void delete_RATcLcLDBLink(void *p);
   static void deleteArray_RATcLcLDBLink(void *p);
   static void destruct_RATcLcLDBLink(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DBLink*)
   {
      ::RAT::DBLink *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::DBLink));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DBLink", "RAT/DBLink.hh", 50,
                  typeid(::RAT::DBLink), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLDBLink_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DBLink) );
      instance.SetDelete(&delete_RATcLcLDBLink);
      instance.SetDeleteArray(&deleteArray_RATcLcLDBLink);
      instance.SetDestructor(&destruct_RATcLcLDBLink);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DBLink*)
   {
      return GenerateInitInstanceLocal((::RAT::DBLink*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DBLink*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLDBLink_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::DBLink*)0x0)->GetClass();
      RATcLcLDBLink_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLDBLink_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *simple_ptr_nocopylERATcLcLDBLinkgR_Dictionary();
   static void simple_ptr_nocopylERATcLcLDBLinkgR_TClassManip(TClass*);
   static void *new_simple_ptr_nocopylERATcLcLDBLinkgR(void *p = 0);
   static void *newArray_simple_ptr_nocopylERATcLcLDBLinkgR(Long_t size, void *p);
   static void delete_simple_ptr_nocopylERATcLcLDBLinkgR(void *p);
   static void deleteArray_simple_ptr_nocopylERATcLcLDBLinkgR(void *p);
   static void destruct_simple_ptr_nocopylERATcLcLDBLinkgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::simple_ptr_nocopy<RAT::DBLink>*)
   {
      ::simple_ptr_nocopy<RAT::DBLink> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::simple_ptr_nocopy<RAT::DBLink>));
      static ::ROOT::TGenericClassInfo 
         instance("simple_ptr_nocopy<RAT::DBLink>", "RAT/smart_ptr.hpp", 353,
                  typeid(::simple_ptr_nocopy<RAT::DBLink>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &simple_ptr_nocopylERATcLcLDBLinkgR_Dictionary, isa_proxy, 4,
                  sizeof(::simple_ptr_nocopy<RAT::DBLink>) );
      instance.SetNew(&new_simple_ptr_nocopylERATcLcLDBLinkgR);
      instance.SetNewArray(&newArray_simple_ptr_nocopylERATcLcLDBLinkgR);
      instance.SetDelete(&delete_simple_ptr_nocopylERATcLcLDBLinkgR);
      instance.SetDeleteArray(&deleteArray_simple_ptr_nocopylERATcLcLDBLinkgR);
      instance.SetDestructor(&destruct_simple_ptr_nocopylERATcLcLDBLinkgR);

      ::ROOT::AddClassAlternate("simple_ptr_nocopy<RAT::DBLink>","RAT::DBLinkPtr");
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::simple_ptr_nocopy<RAT::DBLink>*)
   {
      return GenerateInitInstanceLocal((::simple_ptr_nocopy<RAT::DBLink>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::simple_ptr_nocopy<RAT::DBLink>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *simple_ptr_nocopylERATcLcLDBLinkgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::simple_ptr_nocopy<RAT::DBLink>*)0x0)->GetClass();
      simple_ptr_nocopylERATcLcLDBLinkgR_TClassManip(theClass);
   return theClass;
   }

   static void simple_ptr_nocopylERATcLcLDBLinkgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLDB_Dictionary();
   static void RATcLcLDB_TClassManip(TClass*);
   static void *new_RATcLcLDB(void *p = 0);
   static void *newArray_RATcLcLDB(Long_t size, void *p);
   static void delete_RATcLcLDB(void *p);
   static void deleteArray_RATcLcLDB(void *p);
   static void destruct_RATcLcLDB(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DB*)
   {
      ::RAT::DB *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::DB));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DB", "RAT/DB.hh", 125,
                  typeid(::RAT::DB), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLDB_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DB) );
      instance.SetNew(&new_RATcLcLDB);
      instance.SetNewArray(&newArray_RATcLcLDB);
      instance.SetDelete(&delete_RATcLcLDB);
      instance.SetDeleteArray(&deleteArray_RATcLcLDB);
      instance.SetDestructor(&destruct_RATcLcLDB);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DB*)
   {
      return GenerateInitInstanceLocal((::RAT::DB*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DB*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLDB_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::DB*)0x0)->GetClass();
      RATcLcLDB_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLDB_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLDBTextLoader_Dictionary();
   static void RATcLcLDBTextLoader_TClassManip(TClass*);
   static void *new_RATcLcLDBTextLoader(void *p = 0);
   static void *newArray_RATcLcLDBTextLoader(Long_t size, void *p);
   static void delete_RATcLcLDBTextLoader(void *p);
   static void deleteArray_RATcLcLDBTextLoader(void *p);
   static void destruct_RATcLcLDBTextLoader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DBTextLoader*)
   {
      ::RAT::DBTextLoader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::DBTextLoader));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DBTextLoader", "RAT/DBTextLoader.hh", 105,
                  typeid(::RAT::DBTextLoader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLDBTextLoader_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DBTextLoader) );
      instance.SetNew(&new_RATcLcLDBTextLoader);
      instance.SetNewArray(&newArray_RATcLcLDBTextLoader);
      instance.SetDelete(&delete_RATcLcLDBTextLoader);
      instance.SetDeleteArray(&deleteArray_RATcLcLDBTextLoader);
      instance.SetDestructor(&destruct_RATcLcLDBTextLoader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DBTextLoader*)
   {
      return GenerateInitInstanceLocal((::RAT::DBTextLoader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DBTextLoader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLDBTextLoader_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::DBTextLoader*)0x0)->GetClass();
      RATcLcLDBTextLoader_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLDBTextLoader_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RATcLcLDBJsonLoader_Dictionary();
   static void RATcLcLDBJsonLoader_TClassManip(TClass*);
   static void *new_RATcLcLDBJsonLoader(void *p = 0);
   static void *newArray_RATcLcLDBJsonLoader(Long_t size, void *p);
   static void delete_RATcLcLDBJsonLoader(void *p);
   static void deleteArray_RATcLcLDBJsonLoader(void *p);
   static void destruct_RATcLcLDBJsonLoader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::DBJsonLoader*)
   {
      ::RAT::DBJsonLoader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RAT::DBJsonLoader));
      static ::ROOT::TGenericClassInfo 
         instance("RAT::DBJsonLoader", "RAT/DBJsonLoader.hh", 15,
                  typeid(::RAT::DBJsonLoader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &RATcLcLDBJsonLoader_Dictionary, isa_proxy, 4,
                  sizeof(::RAT::DBJsonLoader) );
      instance.SetNew(&new_RATcLcLDBJsonLoader);
      instance.SetNewArray(&newArray_RATcLcLDBJsonLoader);
      instance.SetDelete(&delete_RATcLcLDBJsonLoader);
      instance.SetDeleteArray(&deleteArray_RATcLcLDBJsonLoader);
      instance.SetDestructor(&destruct_RATcLcLDBJsonLoader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::DBJsonLoader*)
   {
      return GenerateInitInstanceLocal((::RAT::DBJsonLoader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::DBJsonLoader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RATcLcLDBJsonLoader_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RAT::DBJsonLoader*)0x0)->GetClass();
      RATcLcLDBJsonLoader_TClassManip(theClass);
   return theClass;
   }

   static void RATcLcLDBJsonLoader_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLObjInt(void *p = 0);
   static void *newArray_RATcLcLObjInt(Long_t size, void *p);
   static void delete_RATcLcLObjInt(void *p);
   static void deleteArray_RATcLcLObjInt(void *p);
   static void destruct_RATcLcLObjInt(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::ObjInt*)
   {
      ::RAT::ObjInt *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::ObjInt >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::ObjInt", ::RAT::ObjInt::Class_Version(), "RAT/ObjInt.hh", 9,
                  typeid(::RAT::ObjInt), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::ObjInt::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::ObjInt) );
      instance.SetNew(&new_RATcLcLObjInt);
      instance.SetNewArray(&newArray_RATcLcLObjInt);
      instance.SetDelete(&delete_RATcLcLObjInt);
      instance.SetDeleteArray(&deleteArray_RATcLcLObjInt);
      instance.SetDestructor(&destruct_RATcLcLObjInt);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::ObjInt*)
   {
      return GenerateInitInstanceLocal((::RAT::ObjInt*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::ObjInt*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RATcLcLObjDbl(void *p = 0);
   static void *newArray_RATcLcLObjDbl(Long_t size, void *p);
   static void delete_RATcLcLObjDbl(void *p);
   static void deleteArray_RATcLcLObjDbl(void *p);
   static void destruct_RATcLcLObjDbl(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RAT::ObjDbl*)
   {
      ::RAT::ObjDbl *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RAT::ObjDbl >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RAT::ObjDbl", ::RAT::ObjDbl::Class_Version(), "RAT/ObjDbl.hh", 9,
                  typeid(::RAT::ObjDbl), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RAT::ObjDbl::Dictionary, isa_proxy, 4,
                  sizeof(::RAT::ObjDbl) );
      instance.SetNew(&new_RATcLcLObjDbl);
      instance.SetNewArray(&newArray_RATcLcLObjDbl);
      instance.SetDelete(&delete_RATcLcLObjDbl);
      instance.SetDeleteArray(&deleteArray_RATcLcLObjDbl);
      instance.SetDestructor(&destruct_RATcLcLObjDbl);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RAT::ObjDbl*)
   {
      return GenerateInitInstanceLocal((::RAT::ObjDbl*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RAT::ObjDbl*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr MCParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCParticle::Class_Name()
{
   return "RAT::DS::MCParticle";
}

//______________________________________________________________________________
const char *MCParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCParticle*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr MCTrackStep::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCTrackStep::Class_Name()
{
   return "RAT::DS::MCTrackStep";
}

//______________________________________________________________________________
const char *MCTrackStep::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrackStep*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCTrackStep::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrackStep*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCTrackStep::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrackStep*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCTrackStep::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrackStep*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr MCTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCTrack::Class_Name()
{
   return "RAT::DS::MCTrack";
}

//______________________________________________________________________________
const char *MCTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCTrack*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr MCSummary::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCSummary::Class_Name()
{
   return "RAT::DS::MCSummary";
}

//______________________________________________________________________________
const char *MCSummary::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCSummary*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCSummary::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCSummary*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCSummary::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCSummary*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCSummary::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCSummary*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr MCPhoton::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCPhoton::Class_Name()
{
   return "RAT::DS::MCPhoton";
}

//______________________________________________________________________________
const char *MCPhoton::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPhoton*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCPhoton::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPhoton*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCPhoton::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPhoton*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCPhoton::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPhoton*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr MCPMT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCPMT::Class_Name()
{
   return "RAT::DS::MCPMT";
}

//______________________________________________________________________________
const char *MCPMT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPMT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCPMT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPMT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCPMT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPMT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCPMT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MCPMT*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr MC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MC::Class_Name()
{
   return "RAT::DS::MC";
}

//______________________________________________________________________________
const char *MC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::MC*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr PMT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PMT::Class_Name()
{
   return "RAT::DS::PMT";
}

//______________________________________________________________________________
const char *PMT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PMT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PMT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PMT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMT*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr PosFit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PosFit::Class_Name()
{
   return "RAT::DS::PosFit";
}

//______________________________________________________________________________
const char *PosFit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PosFit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PosFit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PosFit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PosFit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PosFit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PosFit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PosFit*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr Centroid::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Centroid::Class_Name()
{
   return "RAT::DS::Centroid";
}

//______________________________________________________________________________
const char *Centroid::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Centroid*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Centroid::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Centroid*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Centroid::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Centroid*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Centroid::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Centroid*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr BonsaiFit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *BonsaiFit::Class_Name()
{
   return "RAT::DS::BonsaiFit";
}

//______________________________________________________________________________
const char *BonsaiFit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::BonsaiFit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int BonsaiFit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::BonsaiFit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *BonsaiFit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::BonsaiFit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *BonsaiFit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::BonsaiFit*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr PathFit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PathFit::Class_Name()
{
   return "RAT::DS::PathFit";
}

//______________________________________________________________________________
const char *PathFit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PathFit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PathFit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PathFit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PathFit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PathFit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PathFit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PathFit*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr EV::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EV::Class_Name()
{
   return "RAT::DS::EV";
}

//______________________________________________________________________________
const char *EV::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::EV*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EV::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::EV*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EV::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::EV*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EV::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::EV*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr Calib::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Calib::Class_Name()
{
   return "RAT::DS::Calib";
}

//______________________________________________________________________________
const char *Calib::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Calib*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Calib::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Calib*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Calib::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Calib*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Calib::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Calib*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr Root::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Root::Class_Name()
{
   return "RAT::DS::Root";
}

//______________________________________________________________________________
const char *Root::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Root*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Root::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Root*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Root::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Root*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Root::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Root*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr PMTInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PMTInfo::Class_Name()
{
   return "RAT::DS::PMTInfo";
}

//______________________________________________________________________________
const char *PMTInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMTInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PMTInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMTInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PMTInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMTInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PMTInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::PMTInfo*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr Run::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Run::Class_Name()
{
   return "RAT::DS::Run";
}

//______________________________________________________________________________
const char *Run::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Run*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Run::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Run*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Run::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Run*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Run::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::Run*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
   namespace DS {
//______________________________________________________________________________
atomic_TClass_ptr RunStore::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RunStore::Class_Name()
{
   return "RAT::DS::RunStore";
}

//______________________________________________________________________________
const char *RunStore::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::RunStore*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RunStore::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::RunStore*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RunStore::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::RunStore*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RunStore::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DS::RunStore*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace RAT {
//______________________________________________________________________________
atomic_TClass_ptr DSReader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DSReader::Class_Name()
{
   return "RAT::DSReader";
}

//______________________________________________________________________________
const char *DSReader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSReader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DSReader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSReader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DSReader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSReader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DSReader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSReader*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT
namespace RAT {
//______________________________________________________________________________
atomic_TClass_ptr DSWriter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DSWriter::Class_Name()
{
   return "RAT::DSWriter";
}

//______________________________________________________________________________
const char *DSWriter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSWriter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DSWriter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSWriter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DSWriter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSWriter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DSWriter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::DSWriter*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT
namespace RAT {
//______________________________________________________________________________
atomic_TClass_ptr TrackNode::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TrackNode::Class_Name()
{
   return "RAT::TrackNode";
}

//______________________________________________________________________________
const char *TrackNode::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::TrackNode*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TrackNode::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::TrackNode*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TrackNode::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::TrackNode*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TrackNode::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::TrackNode*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT
namespace RAT {
//______________________________________________________________________________
atomic_TClass_ptr ObjInt::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ObjInt::Class_Name()
{
   return "RAT::ObjInt";
}

//______________________________________________________________________________
const char *ObjInt::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjInt*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ObjInt::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjInt*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ObjInt::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjInt*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ObjInt::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjInt*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT
namespace RAT {
//______________________________________________________________________________
atomic_TClass_ptr ObjDbl::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ObjDbl::Class_Name()
{
   return "RAT::ObjDbl";
}

//______________________________________________________________________________
const char *ObjDbl::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjDbl*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ObjDbl::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjDbl*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ObjDbl::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjDbl*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ObjDbl::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RAT::ObjDbl*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace RAT
namespace ROOT {
   // Wrappers around operator new
   static void *new_pairlEstringcOstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,string> : new pair<string,string>;
   }
   static void *newArray_pairlEstringcOstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) pair<string,string>[nElements] : new pair<string,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairlEstringcOstringgR(void *p) {
      delete ((pair<string,string>*)p);
   }
   static void deleteArray_pairlEstringcOstringgR(void *p) {
      delete [] ((pair<string,string>*)p);
   }
   static void destruct_pairlEstringcOstringgR(void *p) {
      typedef pair<string,string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class pair<string,string>

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void MCParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::MCParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::MCParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::MCParticle::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLMCParticle(void *p) {
      return  p ? new(p) ::RAT::DS::MCParticle : new ::RAT::DS::MCParticle;
   }
   static void *newArray_RATcLcLDScLcLMCParticle(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::MCParticle[nElements] : new ::RAT::DS::MCParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLMCParticle(void *p) {
      delete ((::RAT::DS::MCParticle*)p);
   }
   static void deleteArray_RATcLcLDScLcLMCParticle(void *p) {
      delete [] ((::RAT::DS::MCParticle*)p);
   }
   static void destruct_RATcLcLDScLcLMCParticle(void *p) {
      typedef ::RAT::DS::MCParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::MCParticle

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void MCTrackStep::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::MCTrackStep.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::MCTrackStep::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::MCTrackStep::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLMCTrackStep(void *p) {
      return  p ? new(p) ::RAT::DS::MCTrackStep : new ::RAT::DS::MCTrackStep;
   }
   static void *newArray_RATcLcLDScLcLMCTrackStep(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::MCTrackStep[nElements] : new ::RAT::DS::MCTrackStep[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLMCTrackStep(void *p) {
      delete ((::RAT::DS::MCTrackStep*)p);
   }
   static void deleteArray_RATcLcLDScLcLMCTrackStep(void *p) {
      delete [] ((::RAT::DS::MCTrackStep*)p);
   }
   static void destruct_RATcLcLDScLcLMCTrackStep(void *p) {
      typedef ::RAT::DS::MCTrackStep current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::MCTrackStep

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void MCTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::MCTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::MCTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::MCTrack::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLMCTrack(void *p) {
      return  p ? new(p) ::RAT::DS::MCTrack : new ::RAT::DS::MCTrack;
   }
   static void *newArray_RATcLcLDScLcLMCTrack(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::MCTrack[nElements] : new ::RAT::DS::MCTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLMCTrack(void *p) {
      delete ((::RAT::DS::MCTrack*)p);
   }
   static void deleteArray_RATcLcLDScLcLMCTrack(void *p) {
      delete [] ((::RAT::DS::MCTrack*)p);
   }
   static void destruct_RATcLcLDScLcLMCTrack(void *p) {
      typedef ::RAT::DS::MCTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::MCTrack

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void MCSummary::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::MCSummary.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::MCSummary::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::MCSummary::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLMCSummary(void *p) {
      return  p ? new(p) ::RAT::DS::MCSummary : new ::RAT::DS::MCSummary;
   }
   static void *newArray_RATcLcLDScLcLMCSummary(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::MCSummary[nElements] : new ::RAT::DS::MCSummary[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLMCSummary(void *p) {
      delete ((::RAT::DS::MCSummary*)p);
   }
   static void deleteArray_RATcLcLDScLcLMCSummary(void *p) {
      delete [] ((::RAT::DS::MCSummary*)p);
   }
   static void destruct_RATcLcLDScLcLMCSummary(void *p) {
      typedef ::RAT::DS::MCSummary current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::MCSummary

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void MCPhoton::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::MCPhoton.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::MCPhoton::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::MCPhoton::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLMCPhoton(void *p) {
      return  p ? new(p) ::RAT::DS::MCPhoton : new ::RAT::DS::MCPhoton;
   }
   static void *newArray_RATcLcLDScLcLMCPhoton(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::MCPhoton[nElements] : new ::RAT::DS::MCPhoton[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLMCPhoton(void *p) {
      delete ((::RAT::DS::MCPhoton*)p);
   }
   static void deleteArray_RATcLcLDScLcLMCPhoton(void *p) {
      delete [] ((::RAT::DS::MCPhoton*)p);
   }
   static void destruct_RATcLcLDScLcLMCPhoton(void *p) {
      typedef ::RAT::DS::MCPhoton current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::MCPhoton

namespace ROOT {
   // Wrappers around operator new
   static void *new_jsoncLcLValue(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::json::Value : new ::json::Value;
   }
   static void *newArray_jsoncLcLValue(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::json::Value[nElements] : new ::json::Value[nElements];
   }
   // Wrapper around operator delete
   static void delete_jsoncLcLValue(void *p) {
      delete ((::json::Value*)p);
   }
   static void deleteArray_jsoncLcLValue(void *p) {
      delete [] ((::json::Value*)p);
   }
   static void destruct_jsoncLcLValue(void *p) {
      typedef ::json::Value current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::json::Value

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLLog(void *p) {
      delete ((::RAT::Log*)p);
   }
   static void deleteArray_RATcLcLLog(void *p) {
      delete [] ((::RAT::Log*)p);
   }
   static void destruct_RATcLcLLog(void *p) {
      typedef ::RAT::Log current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::Log

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void MCPMT::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::MCPMT.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::MCPMT::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::MCPMT::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLMCPMT(void *p) {
      return  p ? new(p) ::RAT::DS::MCPMT : new ::RAT::DS::MCPMT;
   }
   static void *newArray_RATcLcLDScLcLMCPMT(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::MCPMT[nElements] : new ::RAT::DS::MCPMT[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLMCPMT(void *p) {
      delete ((::RAT::DS::MCPMT*)p);
   }
   static void deleteArray_RATcLcLDScLcLMCPMT(void *p) {
      delete [] ((::RAT::DS::MCPMT*)p);
   }
   static void destruct_RATcLcLDScLcLMCPMT(void *p) {
      typedef ::RAT::DS::MCPMT current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::MCPMT

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void MC::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::MC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::MC::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::MC::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLMC(void *p) {
      return  p ? new(p) ::RAT::DS::MC : new ::RAT::DS::MC;
   }
   static void *newArray_RATcLcLDScLcLMC(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::MC[nElements] : new ::RAT::DS::MC[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLMC(void *p) {
      delete ((::RAT::DS::MC*)p);
   }
   static void deleteArray_RATcLcLDScLcLMC(void *p) {
      delete [] ((::RAT::DS::MC*)p);
   }
   static void destruct_RATcLcLDScLcLMC(void *p) {
      typedef ::RAT::DS::MC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::MC

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void PMT::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::PMT.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::PMT::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::PMT::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLPMT(void *p) {
      return  p ? new(p) ::RAT::DS::PMT : new ::RAT::DS::PMT;
   }
   static void *newArray_RATcLcLDScLcLPMT(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::PMT[nElements] : new ::RAT::DS::PMT[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLPMT(void *p) {
      delete ((::RAT::DS::PMT*)p);
   }
   static void deleteArray_RATcLcLDScLcLPMT(void *p) {
      delete [] ((::RAT::DS::PMT*)p);
   }
   static void destruct_RATcLcLDScLcLPMT(void *p) {
      typedef ::RAT::DS::PMT current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::PMT

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void PosFit::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::PosFit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::PosFit::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::PosFit::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLPosFit(void *p) {
      delete ((::RAT::DS::PosFit*)p);
   }
   static void deleteArray_RATcLcLDScLcLPosFit(void *p) {
      delete [] ((::RAT::DS::PosFit*)p);
   }
   static void destruct_RATcLcLDScLcLPosFit(void *p) {
      typedef ::RAT::DS::PosFit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::PosFit

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void Centroid::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::Centroid.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::Centroid::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::Centroid::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLCentroid(void *p) {
      return  p ? new(p) ::RAT::DS::Centroid : new ::RAT::DS::Centroid;
   }
   static void *newArray_RATcLcLDScLcLCentroid(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::Centroid[nElements] : new ::RAT::DS::Centroid[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLCentroid(void *p) {
      delete ((::RAT::DS::Centroid*)p);
   }
   static void deleteArray_RATcLcLDScLcLCentroid(void *p) {
      delete [] ((::RAT::DS::Centroid*)p);
   }
   static void destruct_RATcLcLDScLcLCentroid(void *p) {
      typedef ::RAT::DS::Centroid current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::Centroid

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void BonsaiFit::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::BonsaiFit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::BonsaiFit::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::BonsaiFit::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLBonsaiFit(void *p) {
      return  p ? new(p) ::RAT::DS::BonsaiFit : new ::RAT::DS::BonsaiFit;
   }
   static void *newArray_RATcLcLDScLcLBonsaiFit(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::BonsaiFit[nElements] : new ::RAT::DS::BonsaiFit[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLBonsaiFit(void *p) {
      delete ((::RAT::DS::BonsaiFit*)p);
   }
   static void deleteArray_RATcLcLDScLcLBonsaiFit(void *p) {
      delete [] ((::RAT::DS::BonsaiFit*)p);
   }
   static void destruct_RATcLcLDScLcLBonsaiFit(void *p) {
      typedef ::RAT::DS::BonsaiFit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::BonsaiFit

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void PathFit::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::PathFit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::PathFit::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::PathFit::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLPathFit(void *p) {
      return  p ? new(p) ::RAT::DS::PathFit : new ::RAT::DS::PathFit;
   }
   static void *newArray_RATcLcLDScLcLPathFit(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::PathFit[nElements] : new ::RAT::DS::PathFit[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLPathFit(void *p) {
      delete ((::RAT::DS::PathFit*)p);
   }
   static void deleteArray_RATcLcLDScLcLPathFit(void *p) {
      delete [] ((::RAT::DS::PathFit*)p);
   }
   static void destruct_RATcLcLDScLcLPathFit(void *p) {
      typedef ::RAT::DS::PathFit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::PathFit

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void EV::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::EV.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::EV::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::EV::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLEV(void *p) {
      return  p ? new(p) ::RAT::DS::EV : new ::RAT::DS::EV;
   }
   static void *newArray_RATcLcLDScLcLEV(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::EV[nElements] : new ::RAT::DS::EV[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLEV(void *p) {
      delete ((::RAT::DS::EV*)p);
   }
   static void deleteArray_RATcLcLDScLcLEV(void *p) {
      delete [] ((::RAT::DS::EV*)p);
   }
   static void destruct_RATcLcLDScLcLEV(void *p) {
      typedef ::RAT::DS::EV current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::EV

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void Calib::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::Calib.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::Calib::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::Calib::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLCalib(void *p) {
      return  p ? new(p) ::RAT::DS::Calib : new ::RAT::DS::Calib;
   }
   static void *newArray_RATcLcLDScLcLCalib(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::Calib[nElements] : new ::RAT::DS::Calib[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLCalib(void *p) {
      delete ((::RAT::DS::Calib*)p);
   }
   static void deleteArray_RATcLcLDScLcLCalib(void *p) {
      delete [] ((::RAT::DS::Calib*)p);
   }
   static void destruct_RATcLcLDScLcLCalib(void *p) {
      typedef ::RAT::DS::Calib current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::Calib

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void Root::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::Root.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::Root::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::Root::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLRoot(void *p) {
      return  p ? new(p) ::RAT::DS::Root : new ::RAT::DS::Root;
   }
   static void *newArray_RATcLcLDScLcLRoot(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::Root[nElements] : new ::RAT::DS::Root[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLRoot(void *p) {
      delete ((::RAT::DS::Root*)p);
   }
   static void deleteArray_RATcLcLDScLcLRoot(void *p) {
      delete [] ((::RAT::DS::Root*)p);
   }
   static void destruct_RATcLcLDScLcLRoot(void *p) {
      typedef ::RAT::DS::Root current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::Root

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void PMTInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::PMTInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::PMTInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::PMTInfo::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLPMTInfo(void *p) {
      return  p ? new(p) ::RAT::DS::PMTInfo : new ::RAT::DS::PMTInfo;
   }
   static void *newArray_RATcLcLDScLcLPMTInfo(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::PMTInfo[nElements] : new ::RAT::DS::PMTInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLPMTInfo(void *p) {
      delete ((::RAT::DS::PMTInfo*)p);
   }
   static void deleteArray_RATcLcLDScLcLPMTInfo(void *p) {
      delete [] ((::RAT::DS::PMTInfo*)p);
   }
   static void destruct_RATcLcLDScLcLPMTInfo(void *p) {
      typedef ::RAT::DS::PMTInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::PMTInfo

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void Run::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::Run.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::Run::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::Run::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLRun(void *p) {
      return  p ? new(p) ::RAT::DS::Run : new ::RAT::DS::Run;
   }
   static void *newArray_RATcLcLDScLcLRun(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::DS::Run[nElements] : new ::RAT::DS::Run[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLRun(void *p) {
      delete ((::RAT::DS::Run*)p);
   }
   static void deleteArray_RATcLcLDScLcLRun(void *p) {
      delete [] ((::RAT::DS::Run*)p);
   }
   static void destruct_RATcLcLDScLcLRun(void *p) {
      typedef ::RAT::DS::Run current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::Run

namespace RAT {
   namespace DS {
//______________________________________________________________________________
void RunStore::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DS::RunStore.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::DS::RunStore::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::DS::RunStore::Class(),this);
   }
}

} // namespace RAT::DS
} // namespace RAT::DS
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDScLcLRunStore(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DS::RunStore : new ::RAT::DS::RunStore;
   }
   static void *newArray_RATcLcLDScLcLRunStore(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DS::RunStore[nElements] : new ::RAT::DS::RunStore[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDScLcLRunStore(void *p) {
      delete ((::RAT::DS::RunStore*)p);
   }
   static void deleteArray_RATcLcLDScLcLRunStore(void *p) {
      delete [] ((::RAT::DS::RunStore*)p);
   }
   static void destruct_RATcLcLDScLcLRunStore(void *p) {
      typedef ::RAT::DS::RunStore current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DS::RunStore

namespace RAT {
//______________________________________________________________________________
void DSReader::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DSReader.

   TObject::Streamer(R__b);
}

} // namespace RAT
namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLDSReader(void *p) {
      delete ((::RAT::DSReader*)p);
   }
   static void deleteArray_RATcLcLDSReader(void *p) {
      delete [] ((::RAT::DSReader*)p);
   }
   static void destruct_RATcLcLDSReader(void *p) {
      typedef ::RAT::DSReader current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RATcLcLDSReader(TBuffer &buf, void *obj) {
      ((::RAT::DSReader*)obj)->::RAT::DSReader::Streamer(buf);
   }
} // end of namespace ROOT for class ::RAT::DSReader

namespace RAT {
//______________________________________________________________________________
void DSWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::DSWriter.

   TObject::Streamer(R__b);
}

} // namespace RAT
namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLDSWriter(void *p) {
      delete ((::RAT::DSWriter*)p);
   }
   static void deleteArray_RATcLcLDSWriter(void *p) {
      delete [] ((::RAT::DSWriter*)p);
   }
   static void destruct_RATcLcLDSWriter(void *p) {
      typedef ::RAT::DSWriter current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RATcLcLDSWriter(TBuffer &buf, void *obj) {
      ((::RAT::DSWriter*)obj)->::RAT::DSWriter::Streamer(buf);
   }
} // end of namespace ROOT for class ::RAT::DSWriter

namespace RAT {
//______________________________________________________________________________
void TrackNode::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::TrackNode.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::RAT::TrackNode thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef RAT::DS::MCTrackStep baseClass0;
      baseClass0::Streamer(R__b);
      {
         vector<RAT::TrackNode*> &R__stl =  child;
         R__stl.clear();
         TClass *R__tcl1 = TBuffer::GetClass(typeid(class RAT::TrackNode *));
         if (R__tcl1==0) {
            Error("child streamer","Missing the TClass object for class RAT::TrackNode *!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            RAT::TrackNode* R__t;
            R__t = (RAT::TrackNode*)R__b.ReadObjectAny(R__tcl1);
            R__stl.push_back(R__t);
         }
      }
      R__b >> trackStart;
      R__b >> trackID;
      R__b >> stepID;
      R__b >> pdgcode;
      { TString R__str; R__str.Streamer(R__b); particleName = R__str.Data(); }
      R__b >> prev;
      R__b >> next;
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef RAT::DS::MCTrackStep baseClass0;
      baseClass0::Streamer(R__b);
      {
         vector<RAT::TrackNode*> &R__stl =  child;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<RAT::TrackNode*>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b << trackStart;
      R__b << trackID;
      R__b << stepID;
      R__b << pdgcode;
      { TString R__str = particleName.c_str(); R__str.Streamer(R__b);}
      R__b << prev;
      R__b << next;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} // namespace RAT
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLTrackNode(void *p) {
      return  p ? new(p) ::RAT::TrackNode : new ::RAT::TrackNode;
   }
   static void *newArray_RATcLcLTrackNode(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::TrackNode[nElements] : new ::RAT::TrackNode[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLTrackNode(void *p) {
      delete ((::RAT::TrackNode*)p);
   }
   static void deleteArray_RATcLcLTrackNode(void *p) {
      delete [] ((::RAT::TrackNode*)p);
   }
   static void destruct_RATcLcLTrackNode(void *p) {
      typedef ::RAT::TrackNode current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RATcLcLTrackNode(TBuffer &buf, void *obj) {
      ((::RAT::TrackNode*)obj)->::RAT::TrackNode::Streamer(buf);
   }
} // end of namespace ROOT for class ::RAT::TrackNode

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLTrackCursor(void *p) {
      delete ((::RAT::TrackCursor*)p);
   }
   static void deleteArray_RATcLcLTrackCursor(void *p) {
      delete [] ((::RAT::TrackCursor*)p);
   }
   static void destruct_RATcLcLTrackCursor(void *p) {
      typedef ::RAT::TrackCursor current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::TrackCursor

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLTrackNav(void *p) {
      delete ((::RAT::TrackNav*)p);
   }
   static void deleteArray_RATcLcLTrackNav(void *p) {
      delete [] ((::RAT::TrackNav*)p);
   }
   static void destruct_RATcLcLTrackNav(void *p) {
      typedef ::RAT::TrackNav current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::TrackNav

namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLHTTPDownloader(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::HTTPDownloader : new ::RAT::HTTPDownloader;
   }
   static void *newArray_RATcLcLHTTPDownloader(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::HTTPDownloader[nElements] : new ::RAT::HTTPDownloader[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLHTTPDownloader(void *p) {
      delete ((::RAT::HTTPDownloader*)p);
   }
   static void deleteArray_RATcLcLHTTPDownloader(void *p) {
      delete [] ((::RAT::HTTPDownloader*)p);
   }
   static void destruct_RATcLcLHTTPDownloader(void *p) {
      typedef ::RAT::HTTPDownloader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::HTTPDownloader

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLDBFieldCallback(void *p) {
      delete ((::RAT::DBFieldCallback*)p);
   }
   static void deleteArray_RATcLcLDBFieldCallback(void *p) {
      delete [] ((::RAT::DBFieldCallback*)p);
   }
   static void destruct_RATcLcLDBFieldCallback(void *p) {
      typedef ::RAT::DBFieldCallback current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DBFieldCallback

namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDBTable(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DBTable : new ::RAT::DBTable;
   }
   static void *newArray_RATcLcLDBTable(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DBTable[nElements] : new ::RAT::DBTable[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDBTable(void *p) {
      delete ((::RAT::DBTable*)p);
   }
   static void deleteArray_RATcLcLDBTable(void *p) {
      delete [] ((::RAT::DBTable*)p);
   }
   static void destruct_RATcLcLDBTable(void *p) {
      typedef ::RAT::DBTable current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DBTable

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RATcLcLDBLink(void *p) {
      delete ((::RAT::DBLink*)p);
   }
   static void deleteArray_RATcLcLDBLink(void *p) {
      delete [] ((::RAT::DBLink*)p);
   }
   static void destruct_RATcLcLDBLink(void *p) {
      typedef ::RAT::DBLink current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DBLink

namespace ROOT {
   // Wrappers around operator new
   static void *new_simple_ptr_nocopylERATcLcLDBLinkgR(void *p) {
      return  p ? new(p) ::simple_ptr_nocopy<RAT::DBLink> : new ::simple_ptr_nocopy<RAT::DBLink>;
   }
   static void *newArray_simple_ptr_nocopylERATcLcLDBLinkgR(Long_t nElements, void *p) {
      return p ? new(p) ::simple_ptr_nocopy<RAT::DBLink>[nElements] : new ::simple_ptr_nocopy<RAT::DBLink>[nElements];
   }
   // Wrapper around operator delete
   static void delete_simple_ptr_nocopylERATcLcLDBLinkgR(void *p) {
      delete ((::simple_ptr_nocopy<RAT::DBLink>*)p);
   }
   static void deleteArray_simple_ptr_nocopylERATcLcLDBLinkgR(void *p) {
      delete [] ((::simple_ptr_nocopy<RAT::DBLink>*)p);
   }
   static void destruct_simple_ptr_nocopylERATcLcLDBLinkgR(void *p) {
      typedef ::simple_ptr_nocopy<RAT::DBLink> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::simple_ptr_nocopy<RAT::DBLink>

namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDB(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DB : new ::RAT::DB;
   }
   static void *newArray_RATcLcLDB(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DB[nElements] : new ::RAT::DB[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDB(void *p) {
      delete ((::RAT::DB*)p);
   }
   static void deleteArray_RATcLcLDB(void *p) {
      delete [] ((::RAT::DB*)p);
   }
   static void destruct_RATcLcLDB(void *p) {
      typedef ::RAT::DB current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DB

namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDBTextLoader(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DBTextLoader : new ::RAT::DBTextLoader;
   }
   static void *newArray_RATcLcLDBTextLoader(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DBTextLoader[nElements] : new ::RAT::DBTextLoader[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDBTextLoader(void *p) {
      delete ((::RAT::DBTextLoader*)p);
   }
   static void deleteArray_RATcLcLDBTextLoader(void *p) {
      delete [] ((::RAT::DBTextLoader*)p);
   }
   static void destruct_RATcLcLDBTextLoader(void *p) {
      typedef ::RAT::DBTextLoader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DBTextLoader

namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLDBJsonLoader(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DBJsonLoader : new ::RAT::DBJsonLoader;
   }
   static void *newArray_RATcLcLDBJsonLoader(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::RAT::DBJsonLoader[nElements] : new ::RAT::DBJsonLoader[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLDBJsonLoader(void *p) {
      delete ((::RAT::DBJsonLoader*)p);
   }
   static void deleteArray_RATcLcLDBJsonLoader(void *p) {
      delete [] ((::RAT::DBJsonLoader*)p);
   }
   static void destruct_RATcLcLDBJsonLoader(void *p) {
      typedef ::RAT::DBJsonLoader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::DBJsonLoader

namespace RAT {
//______________________________________________________________________________
void ObjInt::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::ObjInt.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::ObjInt::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::ObjInt::Class(),this);
   }
}

} // namespace RAT
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLObjInt(void *p) {
      return  p ? new(p) ::RAT::ObjInt : new ::RAT::ObjInt;
   }
   static void *newArray_RATcLcLObjInt(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::ObjInt[nElements] : new ::RAT::ObjInt[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLObjInt(void *p) {
      delete ((::RAT::ObjInt*)p);
   }
   static void deleteArray_RATcLcLObjInt(void *p) {
      delete [] ((::RAT::ObjInt*)p);
   }
   static void destruct_RATcLcLObjInt(void *p) {
      typedef ::RAT::ObjInt current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::ObjInt

namespace RAT {
//______________________________________________________________________________
void ObjDbl::Streamer(TBuffer &R__b)
{
   // Stream an object of class RAT::ObjDbl.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RAT::ObjDbl::Class(),this);
   } else {
      R__b.WriteClassBuffer(RAT::ObjDbl::Class(),this);
   }
}

} // namespace RAT
namespace ROOT {
   // Wrappers around operator new
   static void *new_RATcLcLObjDbl(void *p) {
      return  p ? new(p) ::RAT::ObjDbl : new ::RAT::ObjDbl;
   }
   static void *newArray_RATcLcLObjDbl(Long_t nElements, void *p) {
      return p ? new(p) ::RAT::ObjDbl[nElements] : new ::RAT::ObjDbl[nElements];
   }
   // Wrapper around operator delete
   static void delete_RATcLcLObjDbl(void *p) {
      delete ((::RAT::ObjDbl*)p);
   }
   static void deleteArray_RATcLcLObjDbl(void *p) {
      delete [] ((::RAT::ObjDbl*)p);
   }
   static void destruct_RATcLcLObjDbl(void *p) {
      typedef ::RAT::ObjDbl current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RAT::ObjDbl

namespace ROOT {
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary();
   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p);
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<double> >*)
   {
      vector<vector<double> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<double> >", -2, "vector", 386,
                  typeid(vector<vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEdoublegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<double> >) );
      instance.SetNew(&new_vectorlEvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<double> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<double> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<double> >*)0x0)->GetClass();
      vectorlEvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> > : new vector<vector<double> >;
   }
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> >[nElements] : new vector<vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete ((vector<vector<double> >*)p);
   }
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete [] ((vector<vector<double> >*)p);
   }
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p) {
      typedef vector<vector<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<double> >

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 386,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEpairlEstringcOintgRsPgR_Dictionary();
   static void vectorlEpairlEstringcOintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEstringcOintgRsPgR(void *p = 0);
   static void *newArray_vectorlEpairlEstringcOintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEstringcOintgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEstringcOintgRsPgR(void *p);
   static void destruct_vectorlEpairlEstringcOintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<string,int> >*)
   {
      vector<pair<string,int> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<string,int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<string,int> >", -2, "vector", 386,
                  typeid(vector<pair<string,int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEstringcOintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<pair<string,int> >) );
      instance.SetNew(&new_vectorlEpairlEstringcOintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEstringcOintgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEstringcOintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEstringcOintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEstringcOintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<string,int> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<pair<string,int> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEstringcOintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<pair<string,int> >*)0x0)->GetClass();
      vectorlEpairlEstringcOintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEstringcOintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEstringcOintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<string,int> > : new vector<pair<string,int> >;
   }
   static void *newArray_vectorlEpairlEstringcOintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<string,int> >[nElements] : new vector<pair<string,int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEstringcOintgRsPgR(void *p) {
      delete ((vector<pair<string,int> >*)p);
   }
   static void deleteArray_vectorlEpairlEstringcOintgRsPgR(void *p) {
      delete [] ((vector<pair<string,int> >*)p);
   }
   static void destruct_vectorlEpairlEstringcOintgRsPgR(void *p) {
      typedef vector<pair<string,int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<pair<string,int> >

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 386,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 386,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETVector3gR_Dictionary();
   static void vectorlETVector3gR_TClassManip(TClass*);
   static void *new_vectorlETVector3gR(void *p = 0);
   static void *newArray_vectorlETVector3gR(Long_t size, void *p);
   static void delete_vectorlETVector3gR(void *p);
   static void deleteArray_vectorlETVector3gR(void *p);
   static void destruct_vectorlETVector3gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector3>*)
   {
      vector<TVector3> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector3>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector3>", -2, "vector", 386,
                  typeid(vector<TVector3>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETVector3gR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TVector3>) );
      instance.SetNew(&new_vectorlETVector3gR);
      instance.SetNewArray(&newArray_vectorlETVector3gR);
      instance.SetDelete(&delete_vectorlETVector3gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector3gR);
      instance.SetDestructor(&destruct_vectorlETVector3gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector3> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TVector3>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector3gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TVector3>*)0x0)->GetClass();
      vectorlETVector3gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector3gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector3gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector3> : new vector<TVector3>;
   }
   static void *newArray_vectorlETVector3gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector3>[nElements] : new vector<TVector3>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector3gR(void *p) {
      delete ((vector<TVector3>*)p);
   }
   static void deleteArray_vectorlETVector3gR(void *p) {
      delete [] ((vector<TVector3>*)p);
   }
   static void destruct_vectorlETVector3gR(void *p) {
      typedef vector<TVector3> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TVector3>

namespace ROOT {
   static TClass *vectorlERATcLcLTrackNodemUgR_Dictionary();
   static void vectorlERATcLcLTrackNodemUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLTrackNodemUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLTrackNodemUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLTrackNodemUgR(void *p);
   static void deleteArray_vectorlERATcLcLTrackNodemUgR(void *p);
   static void destruct_vectorlERATcLcLTrackNodemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::TrackNode*>*)
   {
      vector<RAT::TrackNode*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::TrackNode*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::TrackNode*>", -2, "vector", 386,
                  typeid(vector<RAT::TrackNode*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLTrackNodemUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::TrackNode*>) );
      instance.SetNew(&new_vectorlERATcLcLTrackNodemUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLTrackNodemUgR);
      instance.SetDelete(&delete_vectorlERATcLcLTrackNodemUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLTrackNodemUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLTrackNodemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::TrackNode*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::TrackNode*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLTrackNodemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::TrackNode*>*)0x0)->GetClass();
      vectorlERATcLcLTrackNodemUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLTrackNodemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLTrackNodemUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::TrackNode*> : new vector<RAT::TrackNode*>;
   }
   static void *newArray_vectorlERATcLcLTrackNodemUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::TrackNode*>[nElements] : new vector<RAT::TrackNode*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLTrackNodemUgR(void *p) {
      delete ((vector<RAT::TrackNode*>*)p);
   }
   static void deleteArray_vectorlERATcLcLTrackNodemUgR(void *p) {
      delete [] ((vector<RAT::TrackNode*>*)p);
   }
   static void destruct_vectorlERATcLcLTrackNodemUgR(void *p) {
      typedef vector<RAT::TrackNode*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::TrackNode*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLRootgR_Dictionary();
   static void vectorlERATcLcLDScLcLRootgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLRootgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLRootgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLRootgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLRootgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLRootgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::Root>*)
   {
      vector<RAT::DS::Root> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::Root>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::Root>", -2, "vector", 386,
                  typeid(vector<RAT::DS::Root>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLRootgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::Root>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLRootgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLRootgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLRootgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLRootgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLRootgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::Root> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::Root>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLRootgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::Root>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLRootgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLRootgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLRootgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::Root> : new vector<RAT::DS::Root>;
   }
   static void *newArray_vectorlERATcLcLDScLcLRootgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::Root>[nElements] : new vector<RAT::DS::Root>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLRootgR(void *p) {
      delete ((vector<RAT::DS::Root>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLRootgR(void *p) {
      delete [] ((vector<RAT::DS::Root>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLRootgR(void *p) {
      typedef vector<RAT::DS::Root> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::Root>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLPathFitgR_Dictionary();
   static void vectorlERATcLcLDScLcLPathFitgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLPathFitgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLPathFitgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLPathFitgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLPathFitgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLPathFitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::PathFit>*)
   {
      vector<RAT::DS::PathFit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::PathFit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::PathFit>", -2, "vector", 386,
                  typeid(vector<RAT::DS::PathFit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLPathFitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RAT::DS::PathFit>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLPathFitgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLPathFitgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLPathFitgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLPathFitgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLPathFitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::PathFit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::PathFit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLPathFitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::PathFit>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLPathFitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLPathFitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLPathFitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PathFit> : new vector<RAT::DS::PathFit>;
   }
   static void *newArray_vectorlERATcLcLDScLcLPathFitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PathFit>[nElements] : new vector<RAT::DS::PathFit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLPathFitgR(void *p) {
      delete ((vector<RAT::DS::PathFit>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLPathFitgR(void *p) {
      delete [] ((vector<RAT::DS::PathFit>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLPathFitgR(void *p) {
      typedef vector<RAT::DS::PathFit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::PathFit>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLPMTInfogR_Dictionary();
   static void vectorlERATcLcLDScLcLPMTInfogR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLPMTInfogR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLPMTInfogR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLPMTInfogR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLPMTInfogR(void *p);
   static void destruct_vectorlERATcLcLDScLcLPMTInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::PMTInfo>*)
   {
      vector<RAT::DS::PMTInfo> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::PMTInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::PMTInfo>", -2, "vector", 386,
                  typeid(vector<RAT::DS::PMTInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLPMTInfogR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::PMTInfo>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLPMTInfogR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLPMTInfogR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLPMTInfogR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLPMTInfogR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLPMTInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::PMTInfo> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::PMTInfo>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLPMTInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::PMTInfo>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLPMTInfogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLPMTInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLPMTInfogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PMTInfo> : new vector<RAT::DS::PMTInfo>;
   }
   static void *newArray_vectorlERATcLcLDScLcLPMTInfogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PMTInfo>[nElements] : new vector<RAT::DS::PMTInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLPMTInfogR(void *p) {
      delete ((vector<RAT::DS::PMTInfo>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLPMTInfogR(void *p) {
      delete [] ((vector<RAT::DS::PMTInfo>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLPMTInfogR(void *p) {
      typedef vector<RAT::DS::PMTInfo> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::PMTInfo>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLPMTgR_Dictionary();
   static void vectorlERATcLcLDScLcLPMTgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLPMTgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLPMTgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLPMTgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLPMTgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLPMTgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::PMT>*)
   {
      vector<RAT::DS::PMT> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::PMT>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::PMT>", -2, "vector", 386,
                  typeid(vector<RAT::DS::PMT>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLPMTgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::PMT>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLPMTgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLPMTgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLPMTgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLPMTgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLPMTgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::PMT> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::PMT>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLPMTgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::PMT>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLPMTgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLPMTgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLPMTgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PMT> : new vector<RAT::DS::PMT>;
   }
   static void *newArray_vectorlERATcLcLDScLcLPMTgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PMT>[nElements] : new vector<RAT::DS::PMT>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLPMTgR(void *p) {
      delete ((vector<RAT::DS::PMT>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLPMTgR(void *p) {
      delete [] ((vector<RAT::DS::PMT>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLPMTgR(void *p) {
      typedef vector<RAT::DS::PMT> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::PMT>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLPMTmUgR_Dictionary();
   static void vectorlERATcLcLDScLcLPMTmUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLPMTmUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLPMTmUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLPMTmUgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLPMTmUgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLPMTmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::PMT*>*)
   {
      vector<RAT::DS::PMT*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::PMT*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::PMT*>", -2, "vector", 386,
                  typeid(vector<RAT::DS::PMT*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLPMTmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::PMT*>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLPMTmUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLPMTmUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLPMTmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLPMTmUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLPMTmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::PMT*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::PMT*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLPMTmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::PMT*>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLPMTmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLPMTmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLPMTmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PMT*> : new vector<RAT::DS::PMT*>;
   }
   static void *newArray_vectorlERATcLcLDScLcLPMTmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::PMT*>[nElements] : new vector<RAT::DS::PMT*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLPMTmUgR(void *p) {
      delete ((vector<RAT::DS::PMT*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLPMTmUgR(void *p) {
      delete [] ((vector<RAT::DS::PMT*>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLPMTmUgR(void *p) {
      typedef vector<RAT::DS::PMT*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::PMT*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCTrackStepgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCTrackStepgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCTrackStepgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCTrackStepgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCTrackStepgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackStepgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCTrackStepgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCTrackStep>*)
   {
      vector<RAT::DS::MCTrackStep> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCTrackStep>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCTrackStep>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCTrackStep>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCTrackStepgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCTrackStep>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCTrackStepgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCTrackStepgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCTrackStepgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCTrackStepgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCTrackStepgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCTrackStep> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCTrackStep>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCTrackStepgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCTrackStep>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCTrackStepgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCTrackStepgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCTrackStepgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrackStep> : new vector<RAT::DS::MCTrackStep>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCTrackStepgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrackStep>[nElements] : new vector<RAT::DS::MCTrackStep>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCTrackStepgR(void *p) {
      delete ((vector<RAT::DS::MCTrackStep>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackStepgR(void *p) {
      delete [] ((vector<RAT::DS::MCTrackStep>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCTrackStepgR(void *p) {
      typedef vector<RAT::DS::MCTrackStep> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCTrackStep>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCTrackStepmUgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCTrackStepmUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCTrackStepmUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCTrackStep*>*)
   {
      vector<RAT::DS::MCTrackStep*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCTrackStep*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCTrackStep*>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCTrackStep*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCTrackStepmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCTrackStep*>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCTrackStepmUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCTrackStepmUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCTrackStepmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCTrackStepmUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCTrackStepmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCTrackStep*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCTrackStep*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCTrackStepmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCTrackStep*>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCTrackStepmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCTrackStepmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrackStep*> : new vector<RAT::DS::MCTrackStep*>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCTrackStepmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrackStep*>[nElements] : new vector<RAT::DS::MCTrackStep*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p) {
      delete ((vector<RAT::DS::MCTrackStep*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p) {
      delete [] ((vector<RAT::DS::MCTrackStep*>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCTrackStepmUgR(void *p) {
      typedef vector<RAT::DS::MCTrackStep*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCTrackStep*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCTrackgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCTrackgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCTrackgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCTrackgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCTrackgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCTrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCTrack>*)
   {
      vector<RAT::DS::MCTrack> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCTrack>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCTrack>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCTrack>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCTrackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCTrack>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCTrackgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCTrackgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCTrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCTrackgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCTrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCTrack> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCTrack>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCTrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCTrack>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCTrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCTrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCTrackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrack> : new vector<RAT::DS::MCTrack>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCTrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrack>[nElements] : new vector<RAT::DS::MCTrack>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCTrackgR(void *p) {
      delete ((vector<RAT::DS::MCTrack>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackgR(void *p) {
      delete [] ((vector<RAT::DS::MCTrack>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCTrackgR(void *p) {
      typedef vector<RAT::DS::MCTrack> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCTrack>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCTrackmUgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCTrackmUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCTrackmUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCTrackmUgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackmUgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCTrack*>*)
   {
      vector<RAT::DS::MCTrack*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCTrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCTrack*>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCTrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCTrackmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCTrack*>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCTrackmUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCTrackmUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCTrackmUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCTrack*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCTrack*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCTrack*>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCTrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCTrackmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrack*> : new vector<RAT::DS::MCTrack*>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCTrack*>[nElements] : new vector<RAT::DS::MCTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCTrackmUgR(void *p) {
      delete ((vector<RAT::DS::MCTrack*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCTrackmUgR(void *p) {
      delete [] ((vector<RAT::DS::MCTrack*>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCTrackmUgR(void *p) {
      typedef vector<RAT::DS::MCTrack*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCTrack*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCSummarygR_Dictionary();
   static void vectorlERATcLcLDScLcLMCSummarygR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCSummarygR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCSummarygR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCSummarygR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCSummarygR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCSummarygR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCSummary>*)
   {
      vector<RAT::DS::MCSummary> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCSummary>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCSummary>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCSummary>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCSummarygR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RAT::DS::MCSummary>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCSummarygR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCSummarygR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCSummarygR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCSummarygR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCSummarygR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCSummary> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCSummary>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCSummarygR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCSummary>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCSummarygR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCSummarygR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCSummarygR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCSummary> : new vector<RAT::DS::MCSummary>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCSummarygR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCSummary>[nElements] : new vector<RAT::DS::MCSummary>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCSummarygR(void *p) {
      delete ((vector<RAT::DS::MCSummary>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCSummarygR(void *p) {
      delete [] ((vector<RAT::DS::MCSummary>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCSummarygR(void *p) {
      typedef vector<RAT::DS::MCSummary> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCSummary>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCPhotongR_Dictionary();
   static void vectorlERATcLcLDScLcLMCPhotongR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCPhotongR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCPhotongR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCPhotongR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCPhotongR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCPhotongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCPhoton>*)
   {
      vector<RAT::DS::MCPhoton> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCPhoton>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCPhoton>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCPhoton>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCPhotongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCPhoton>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCPhotongR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCPhotongR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCPhotongR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCPhotongR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCPhotongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCPhoton> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCPhoton>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCPhotongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCPhoton>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCPhotongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCPhotongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCPhotongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPhoton> : new vector<RAT::DS::MCPhoton>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCPhotongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPhoton>[nElements] : new vector<RAT::DS::MCPhoton>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCPhotongR(void *p) {
      delete ((vector<RAT::DS::MCPhoton>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCPhotongR(void *p) {
      delete [] ((vector<RAT::DS::MCPhoton>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCPhotongR(void *p) {
      typedef vector<RAT::DS::MCPhoton> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCPhoton>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCPhotonmUgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCPhotonmUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCPhotonmUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCPhoton*>*)
   {
      vector<RAT::DS::MCPhoton*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCPhoton*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCPhoton*>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCPhoton*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCPhotonmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCPhoton*>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCPhotonmUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCPhotonmUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCPhotonmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCPhotonmUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCPhotonmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCPhoton*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCPhoton*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCPhotonmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCPhoton*>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCPhotonmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCPhotonmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPhoton*> : new vector<RAT::DS::MCPhoton*>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCPhotonmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPhoton*>[nElements] : new vector<RAT::DS::MCPhoton*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p) {
      delete ((vector<RAT::DS::MCPhoton*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p) {
      delete [] ((vector<RAT::DS::MCPhoton*>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCPhotonmUgR(void *p) {
      typedef vector<RAT::DS::MCPhoton*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCPhoton*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCParticlegR_Dictionary();
   static void vectorlERATcLcLDScLcLMCParticlegR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCParticlegR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCParticlegR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCParticlegR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCParticlegR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCParticle>*)
   {
      vector<RAT::DS::MCParticle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCParticle>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCParticle>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCParticlegR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCParticlegR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCParticlegR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCParticle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCParticle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCParticle>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCParticle> : new vector<RAT::DS::MCParticle>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCParticle>[nElements] : new vector<RAT::DS::MCParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCParticlegR(void *p) {
      delete ((vector<RAT::DS::MCParticle>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCParticlegR(void *p) {
      delete [] ((vector<RAT::DS::MCParticle>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCParticlegR(void *p) {
      typedef vector<RAT::DS::MCParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCParticle>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCParticlemUgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCParticlemUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCParticlemUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCParticlemUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCParticlemUgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCParticlemUgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCParticlemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCParticle*>*)
   {
      vector<RAT::DS::MCParticle*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCParticle*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCParticle*>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCParticle*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCParticlemUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCParticle*>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCParticlemUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCParticlemUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCParticlemUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCParticlemUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCParticlemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCParticle*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCParticle*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCParticlemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCParticle*>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCParticlemUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCParticlemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCParticlemUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCParticle*> : new vector<RAT::DS::MCParticle*>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCParticlemUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCParticle*>[nElements] : new vector<RAT::DS::MCParticle*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCParticlemUgR(void *p) {
      delete ((vector<RAT::DS::MCParticle*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCParticlemUgR(void *p) {
      delete [] ((vector<RAT::DS::MCParticle*>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCParticlemUgR(void *p) {
      typedef vector<RAT::DS::MCParticle*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCParticle*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCPMTgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCPMTgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCPMTgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCPMTgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCPMTgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCPMTgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCPMTgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCPMT>*)
   {
      vector<RAT::DS::MCPMT> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCPMT>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCPMT>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCPMT>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCPMTgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCPMT>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCPMTgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCPMTgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCPMTgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCPMTgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCPMTgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCPMT> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCPMT>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCPMTgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCPMT>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCPMTgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCPMTgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCPMTgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPMT> : new vector<RAT::DS::MCPMT>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCPMTgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPMT>[nElements] : new vector<RAT::DS::MCPMT>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCPMTgR(void *p) {
      delete ((vector<RAT::DS::MCPMT>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCPMTgR(void *p) {
      delete [] ((vector<RAT::DS::MCPMT>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCPMTgR(void *p) {
      typedef vector<RAT::DS::MCPMT> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCPMT>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCPMTmUgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCPMTmUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCPMTmUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCPMTmUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCPMTmUgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCPMTmUgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCPMTmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MCPMT*>*)
   {
      vector<RAT::DS::MCPMT*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MCPMT*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MCPMT*>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MCPMT*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCPMTmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MCPMT*>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCPMTmUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCPMTmUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCPMTmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCPMTmUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCPMTmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MCPMT*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MCPMT*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCPMTmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MCPMT*>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCPMTmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCPMTmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCPMTmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPMT*> : new vector<RAT::DS::MCPMT*>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCPMTmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MCPMT*>[nElements] : new vector<RAT::DS::MCPMT*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCPMTmUgR(void *p) {
      delete ((vector<RAT::DS::MCPMT*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCPMTmUgR(void *p) {
      delete [] ((vector<RAT::DS::MCPMT*>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCPMTmUgR(void *p) {
      typedef vector<RAT::DS::MCPMT*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MCPMT*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLMCgR_Dictionary();
   static void vectorlERATcLcLDScLcLMCgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLMCgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLMCgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLMCgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLMCgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLMCgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::MC>*)
   {
      vector<RAT::DS::MC> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::MC>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::MC>", -2, "vector", 386,
                  typeid(vector<RAT::DS::MC>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLMCgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::MC>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLMCgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLMCgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLMCgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLMCgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLMCgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::MC> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::MC>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLMCgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::MC>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLMCgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLMCgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLMCgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MC> : new vector<RAT::DS::MC>;
   }
   static void *newArray_vectorlERATcLcLDScLcLMCgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::MC>[nElements] : new vector<RAT::DS::MC>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLMCgR(void *p) {
      delete ((vector<RAT::DS::MC>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLMCgR(void *p) {
      delete [] ((vector<RAT::DS::MC>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLMCgR(void *p) {
      typedef vector<RAT::DS::MC> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::MC>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLEVgR_Dictionary();
   static void vectorlERATcLcLDScLcLEVgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLEVgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLEVgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLEVgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLEVgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLEVgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::EV>*)
   {
      vector<RAT::DS::EV> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::EV>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::EV>", -2, "vector", 386,
                  typeid(vector<RAT::DS::EV>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLEVgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::EV>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLEVgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLEVgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLEVgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLEVgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLEVgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::EV> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::EV>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLEVgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::EV>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLEVgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLEVgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLEVgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::EV> : new vector<RAT::DS::EV>;
   }
   static void *newArray_vectorlERATcLcLDScLcLEVgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::EV>[nElements] : new vector<RAT::DS::EV>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLEVgR(void *p) {
      delete ((vector<RAT::DS::EV>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLEVgR(void *p) {
      delete [] ((vector<RAT::DS::EV>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLEVgR(void *p) {
      typedef vector<RAT::DS::EV> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::EV>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLEVmUgR_Dictionary();
   static void vectorlERATcLcLDScLcLEVmUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLEVmUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLEVmUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLEVmUgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLEVmUgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLEVmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::EV*>*)
   {
      vector<RAT::DS::EV*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::EV*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::EV*>", -2, "vector", 386,
                  typeid(vector<RAT::DS::EV*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLEVmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::EV*>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLEVmUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLEVmUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLEVmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLEVmUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLEVmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::EV*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::EV*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLEVmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::EV*>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLEVmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLEVmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLEVmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::EV*> : new vector<RAT::DS::EV*>;
   }
   static void *newArray_vectorlERATcLcLDScLcLEVmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::EV*>[nElements] : new vector<RAT::DS::EV*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLEVmUgR(void *p) {
      delete ((vector<RAT::DS::EV*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLEVmUgR(void *p) {
      delete [] ((vector<RAT::DS::EV*>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLEVmUgR(void *p) {
      typedef vector<RAT::DS::EV*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::EV*>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLCentroidgR_Dictionary();
   static void vectorlERATcLcLDScLcLCentroidgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLCentroidgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLCentroidgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLCentroidgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLCentroidgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLCentroidgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::Centroid>*)
   {
      vector<RAT::DS::Centroid> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::Centroid>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::Centroid>", -2, "vector", 386,
                  typeid(vector<RAT::DS::Centroid>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLCentroidgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RAT::DS::Centroid>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLCentroidgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLCentroidgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLCentroidgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLCentroidgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLCentroidgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::Centroid> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::Centroid>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLCentroidgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::Centroid>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLCentroidgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLCentroidgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLCentroidgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::Centroid> : new vector<RAT::DS::Centroid>;
   }
   static void *newArray_vectorlERATcLcLDScLcLCentroidgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::Centroid>[nElements] : new vector<RAT::DS::Centroid>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLCentroidgR(void *p) {
      delete ((vector<RAT::DS::Centroid>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLCentroidgR(void *p) {
      delete [] ((vector<RAT::DS::Centroid>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLCentroidgR(void *p) {
      typedef vector<RAT::DS::Centroid> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::Centroid>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLCalibgR_Dictionary();
   static void vectorlERATcLcLDScLcLCalibgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLCalibgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLCalibgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLCalibgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLCalibgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLCalibgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::Calib>*)
   {
      vector<RAT::DS::Calib> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::Calib>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::Calib>", -2, "vector", 386,
                  typeid(vector<RAT::DS::Calib>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLCalibgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DS::Calib>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLCalibgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLCalibgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLCalibgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLCalibgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLCalibgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::Calib> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::Calib>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLCalibgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::Calib>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLCalibgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLCalibgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLCalibgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::Calib> : new vector<RAT::DS::Calib>;
   }
   static void *newArray_vectorlERATcLcLDScLcLCalibgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::Calib>[nElements] : new vector<RAT::DS::Calib>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLCalibgR(void *p) {
      delete ((vector<RAT::DS::Calib>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLCalibgR(void *p) {
      delete [] ((vector<RAT::DS::Calib>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLCalibgR(void *p) {
      typedef vector<RAT::DS::Calib> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::Calib>

namespace ROOT {
   static TClass *vectorlERATcLcLDScLcLBonsaiFitgR_Dictionary();
   static void vectorlERATcLcLDScLcLBonsaiFitgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDScLcLBonsaiFitgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDScLcLBonsaiFitgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDScLcLBonsaiFitgR(void *p);
   static void deleteArray_vectorlERATcLcLDScLcLBonsaiFitgR(void *p);
   static void destruct_vectorlERATcLcLDScLcLBonsaiFitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DS::BonsaiFit>*)
   {
      vector<RAT::DS::BonsaiFit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DS::BonsaiFit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DS::BonsaiFit>", -2, "vector", 386,
                  typeid(vector<RAT::DS::BonsaiFit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDScLcLBonsaiFitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RAT::DS::BonsaiFit>) );
      instance.SetNew(&new_vectorlERATcLcLDScLcLBonsaiFitgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDScLcLBonsaiFitgR);
      instance.SetDelete(&delete_vectorlERATcLcLDScLcLBonsaiFitgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDScLcLBonsaiFitgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDScLcLBonsaiFitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DS::BonsaiFit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DS::BonsaiFit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDScLcLBonsaiFitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DS::BonsaiFit>*)0x0)->GetClass();
      vectorlERATcLcLDScLcLBonsaiFitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDScLcLBonsaiFitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDScLcLBonsaiFitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::BonsaiFit> : new vector<RAT::DS::BonsaiFit>;
   }
   static void *newArray_vectorlERATcLcLDScLcLBonsaiFitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DS::BonsaiFit>[nElements] : new vector<RAT::DS::BonsaiFit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDScLcLBonsaiFitgR(void *p) {
      delete ((vector<RAT::DS::BonsaiFit>*)p);
   }
   static void deleteArray_vectorlERATcLcLDScLcLBonsaiFitgR(void *p) {
      delete [] ((vector<RAT::DS::BonsaiFit>*)p);
   }
   static void destruct_vectorlERATcLcLDScLcLBonsaiFitgR(void *p) {
      typedef vector<RAT::DS::BonsaiFit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DS::BonsaiFit>

namespace ROOT {
   static TClass *vectorlERATcLcLDBTablemUgR_Dictionary();
   static void vectorlERATcLcLDBTablemUgR_TClassManip(TClass*);
   static void *new_vectorlERATcLcLDBTablemUgR(void *p = 0);
   static void *newArray_vectorlERATcLcLDBTablemUgR(Long_t size, void *p);
   static void delete_vectorlERATcLcLDBTablemUgR(void *p);
   static void deleteArray_vectorlERATcLcLDBTablemUgR(void *p);
   static void destruct_vectorlERATcLcLDBTablemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RAT::DBTable*>*)
   {
      vector<RAT::DBTable*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RAT::DBTable*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RAT::DBTable*>", -2, "vector", 386,
                  typeid(vector<RAT::DBTable*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERATcLcLDBTablemUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RAT::DBTable*>) );
      instance.SetNew(&new_vectorlERATcLcLDBTablemUgR);
      instance.SetNewArray(&newArray_vectorlERATcLcLDBTablemUgR);
      instance.SetDelete(&delete_vectorlERATcLcLDBTablemUgR);
      instance.SetDeleteArray(&deleteArray_vectorlERATcLcLDBTablemUgR);
      instance.SetDestructor(&destruct_vectorlERATcLcLDBTablemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RAT::DBTable*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RAT::DBTable*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERATcLcLDBTablemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RAT::DBTable*>*)0x0)->GetClass();
      vectorlERATcLcLDBTablemUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERATcLcLDBTablemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERATcLcLDBTablemUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DBTable*> : new vector<RAT::DBTable*>;
   }
   static void *newArray_vectorlERATcLcLDBTablemUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RAT::DBTable*>[nElements] : new vector<RAT::DBTable*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERATcLcLDBTablemUgR(void *p) {
      delete ((vector<RAT::DBTable*>*)p);
   }
   static void deleteArray_vectorlERATcLcLDBTablemUgR(void *p) {
      delete [] ((vector<RAT::DBTable*>*)p);
   }
   static void destruct_vectorlERATcLcLDBTablemUgR(void *p) {
      typedef vector<RAT::DBTable*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RAT::DBTable*>

namespace ROOT {
   static TClass *setlEstringgR_Dictionary();
   static void setlEstringgR_TClassManip(TClass*);
   static void *new_setlEstringgR(void *p = 0);
   static void *newArray_setlEstringgR(Long_t size, void *p);
   static void delete_setlEstringgR(void *p);
   static void deleteArray_setlEstringgR(void *p);
   static void destruct_setlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<string>*)
   {
      set<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<string>));
      static ::ROOT::TGenericClassInfo 
         instance("set<string>", -2, "set", 94,
                  typeid(set<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(set<string>) );
      instance.SetNew(&new_setlEstringgR);
      instance.SetNewArray(&newArray_setlEstringgR);
      instance.SetDelete(&delete_setlEstringgR);
      instance.SetDeleteArray(&deleteArray_setlEstringgR);
      instance.SetDestructor(&destruct_setlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<string>*)0x0)->GetClass();
      setlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<string> : new set<string>;
   }
   static void *newArray_setlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<string>[nElements] : new set<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEstringgR(void *p) {
      delete ((set<string>*)p);
   }
   static void deleteArray_setlEstringgR(void *p) {
      delete [] ((set<string>*)p);
   }
   static void destruct_setlEstringgR(void *p) {
      typedef set<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<string>

namespace ROOT {
   static TClass *setlERATcLcLDBTableKeygR_Dictionary();
   static void setlERATcLcLDBTableKeygR_TClassManip(TClass*);
   static void *new_setlERATcLcLDBTableKeygR(void *p = 0);
   static void *newArray_setlERATcLcLDBTableKeygR(Long_t size, void *p);
   static void delete_setlERATcLcLDBTableKeygR(void *p);
   static void deleteArray_setlERATcLcLDBTableKeygR(void *p);
   static void destruct_setlERATcLcLDBTableKeygR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<RAT::DBTableKey>*)
   {
      set<RAT::DBTableKey> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<RAT::DBTableKey>));
      static ::ROOT::TGenericClassInfo 
         instance("set<RAT::DBTableKey>", -2, "set", 94,
                  typeid(set<RAT::DBTableKey>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlERATcLcLDBTableKeygR_Dictionary, isa_proxy, 0,
                  sizeof(set<RAT::DBTableKey>) );
      instance.SetNew(&new_setlERATcLcLDBTableKeygR);
      instance.SetNewArray(&newArray_setlERATcLcLDBTableKeygR);
      instance.SetDelete(&delete_setlERATcLcLDBTableKeygR);
      instance.SetDeleteArray(&deleteArray_setlERATcLcLDBTableKeygR);
      instance.SetDestructor(&destruct_setlERATcLcLDBTableKeygR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<RAT::DBTableKey> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<RAT::DBTableKey>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlERATcLcLDBTableKeygR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<RAT::DBTableKey>*)0x0)->GetClass();
      setlERATcLcLDBTableKeygR_TClassManip(theClass);
   return theClass;
   }

   static void setlERATcLcLDBTableKeygR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlERATcLcLDBTableKeygR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<RAT::DBTableKey> : new set<RAT::DBTableKey>;
   }
   static void *newArray_setlERATcLcLDBTableKeygR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<RAT::DBTableKey>[nElements] : new set<RAT::DBTableKey>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlERATcLcLDBTableKeygR(void *p) {
      delete ((set<RAT::DBTableKey>*)p);
   }
   static void deleteArray_setlERATcLcLDBTableKeygR(void *p) {
      delete [] ((set<RAT::DBTableKey>*)p);
   }
   static void destruct_setlERATcLcLDBTableKeygR(void *p) {
      typedef set<RAT::DBTableKey> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<RAT::DBTableKey>

namespace ROOT {
   static TClass *maplEstringcOdoublegR_Dictionary();
   static void maplEstringcOdoublegR_TClassManip(TClass*);
   static void *new_maplEstringcOdoublegR(void *p = 0);
   static void *newArray_maplEstringcOdoublegR(Long_t size, void *p);
   static void delete_maplEstringcOdoublegR(void *p);
   static void deleteArray_maplEstringcOdoublegR(void *p);
   static void destruct_maplEstringcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,double>*)
   {
      map<string,double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,double>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,double>", -2, "map", 100,
                  typeid(map<string,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,double>) );
      instance.SetNew(&new_maplEstringcOdoublegR);
      instance.SetNewArray(&newArray_maplEstringcOdoublegR);
      instance.SetDelete(&delete_maplEstringcOdoublegR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOdoublegR);
      instance.SetDestructor(&destruct_maplEstringcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,double>*)0x0)->GetClass();
      maplEstringcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,double> : new map<string,double>;
   }
   static void *newArray_maplEstringcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,double>[nElements] : new map<string,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOdoublegR(void *p) {
      delete ((map<string,double>*)p);
   }
   static void deleteArray_maplEstringcOdoublegR(void *p) {
      delete [] ((map<string,double>*)p);
   }
   static void destruct_maplEstringcOdoublegR(void *p) {
      typedef map<string,double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,double>

namespace ROOT {
   static TClass *maplEintcOvectorlEintgRsPgR_Dictionary();
   static void maplEintcOvectorlEintgRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOvectorlEintgRsPgR(void *p = 0);
   static void *newArray_maplEintcOvectorlEintgRsPgR(Long_t size, void *p);
   static void delete_maplEintcOvectorlEintgRsPgR(void *p);
   static void deleteArray_maplEintcOvectorlEintgRsPgR(void *p);
   static void destruct_maplEintcOvectorlEintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,vector<int> >*)
   {
      map<int,vector<int> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,vector<int> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,vector<int> >", -2, "map", 100,
                  typeid(map<int,vector<int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOvectorlEintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(map<int,vector<int> >) );
      instance.SetNew(&new_maplEintcOvectorlEintgRsPgR);
      instance.SetNewArray(&newArray_maplEintcOvectorlEintgRsPgR);
      instance.SetDelete(&delete_maplEintcOvectorlEintgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOvectorlEintgRsPgR);
      instance.SetDestructor(&destruct_maplEintcOvectorlEintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,vector<int> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,vector<int> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOvectorlEintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,vector<int> >*)0x0)->GetClass();
      maplEintcOvectorlEintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOvectorlEintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOvectorlEintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,vector<int> > : new map<int,vector<int> >;
   }
   static void *newArray_maplEintcOvectorlEintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,vector<int> >[nElements] : new map<int,vector<int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOvectorlEintgRsPgR(void *p) {
      delete ((map<int,vector<int> >*)p);
   }
   static void deleteArray_maplEintcOvectorlEintgRsPgR(void *p) {
      delete [] ((map<int,vector<int> >*)p);
   }
   static void destruct_maplEintcOvectorlEintgRsPgR(void *p) {
      typedef map<int,vector<int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,vector<int> >

namespace ROOT {
   static TClass *maplEintcOintgR_Dictionary();
   static void maplEintcOintgR_TClassManip(TClass*);
   static void *new_maplEintcOintgR(void *p = 0);
   static void *newArray_maplEintcOintgR(Long_t size, void *p);
   static void delete_maplEintcOintgR(void *p);
   static void deleteArray_maplEintcOintgR(void *p);
   static void destruct_maplEintcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,int>*)
   {
      map<int,int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,int>", -2, "map", 100,
                  typeid(map<int,int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOintgR_Dictionary, isa_proxy, 4,
                  sizeof(map<int,int>) );
      instance.SetNew(&new_maplEintcOintgR);
      instance.SetNewArray(&newArray_maplEintcOintgR);
      instance.SetDelete(&delete_maplEintcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOintgR);
      instance.SetDestructor(&destruct_maplEintcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,int>*)0x0)->GetClass();
      maplEintcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int> : new map<int,int>;
   }
   static void *newArray_maplEintcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,int>[nElements] : new map<int,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOintgR(void *p) {
      delete ((map<int,int>*)p);
   }
   static void deleteArray_maplEintcOintgR(void *p) {
      delete [] ((map<int,int>*)p);
   }
   static void destruct_maplEintcOintgR(void *p) {
      typedef map<int,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,int>

namespace ROOT {
   static TClass *maplEintcORATcLcLDScLcLRunRecordmUgR_Dictionary();
   static void maplEintcORATcLcLDScLcLRunRecordmUgR_TClassManip(TClass*);
   static void *new_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p = 0);
   static void *newArray_maplEintcORATcLcLDScLcLRunRecordmUgR(Long_t size, void *p);
   static void delete_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p);
   static void deleteArray_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p);
   static void destruct_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,RAT::DS::RunRecord*>*)
   {
      map<int,RAT::DS::RunRecord*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,RAT::DS::RunRecord*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,RAT::DS::RunRecord*>", -2, "map", 100,
                  typeid(map<int,RAT::DS::RunRecord*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcORATcLcLDScLcLRunRecordmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,RAT::DS::RunRecord*>) );
      instance.SetNew(&new_maplEintcORATcLcLDScLcLRunRecordmUgR);
      instance.SetNewArray(&newArray_maplEintcORATcLcLDScLcLRunRecordmUgR);
      instance.SetDelete(&delete_maplEintcORATcLcLDScLcLRunRecordmUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcORATcLcLDScLcLRunRecordmUgR);
      instance.SetDestructor(&destruct_maplEintcORATcLcLDScLcLRunRecordmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,RAT::DS::RunRecord*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,RAT::DS::RunRecord*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcORATcLcLDScLcLRunRecordmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,RAT::DS::RunRecord*>*)0x0)->GetClass();
      maplEintcORATcLcLDScLcLRunRecordmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcORATcLcLDScLcLRunRecordmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,RAT::DS::RunRecord*> : new map<int,RAT::DS::RunRecord*>;
   }
   static void *newArray_maplEintcORATcLcLDScLcLRunRecordmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,RAT::DS::RunRecord*>[nElements] : new map<int,RAT::DS::RunRecord*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p) {
      delete ((map<int,RAT::DS::RunRecord*>*)p);
   }
   static void deleteArray_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p) {
      delete [] ((map<int,RAT::DS::RunRecord*>*)p);
   }
   static void destruct_maplEintcORATcLcLDScLcLRunRecordmUgR(void *p) {
      typedef map<int,RAT::DS::RunRecord*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,RAT::DS::RunRecord*>

namespace ROOT {
   static TClass *maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR_Dictionary();
   static void maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR_TClassManip(TClass*);
   static void *new_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p = 0);
   static void *newArray_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(Long_t size, void *p);
   static void delete_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p);
   static void deleteArray_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p);
   static void destruct_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >*)
   {
      map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >", -2, "map", 100,
                  typeid(map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >) );
      instance.SetNew(&new_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR);
      instance.SetNewArray(&newArray_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR);
      instance.SetDelete(&delete_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR);
      instance.SetDeleteArray(&deleteArray_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR);
      instance.SetDestructor(&destruct_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >*)0x0)->GetClass();
      maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> > : new map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >;
   }
   static void *newArray_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >[nElements] : new map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p) {
      delete ((map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >*)p);
   }
   static void deleteArray_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p) {
      delete [] ((map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >*)p);
   }
   static void destruct_maplERATcLcLDBTableKeycOsimple_ptr_nocopylERATcLcLDBTablegRsPgR(void *p) {
      typedef map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<RAT::DBTableKey,simple_ptr_nocopy<RAT::DBTable> >

namespace ROOT {
   static TClass *listlERATcLcLDBLinkmUgR_Dictionary();
   static void listlERATcLcLDBLinkmUgR_TClassManip(TClass*);
   static void *new_listlERATcLcLDBLinkmUgR(void *p = 0);
   static void *newArray_listlERATcLcLDBLinkmUgR(Long_t size, void *p);
   static void delete_listlERATcLcLDBLinkmUgR(void *p);
   static void deleteArray_listlERATcLcLDBLinkmUgR(void *p);
   static void destruct_listlERATcLcLDBLinkmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const list<RAT::DBLink*>*)
   {
      list<RAT::DBLink*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(list<RAT::DBLink*>));
      static ::ROOT::TGenericClassInfo 
         instance("list<RAT::DBLink*>", -2, "list", 552,
                  typeid(list<RAT::DBLink*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &listlERATcLcLDBLinkmUgR_Dictionary, isa_proxy, 0,
                  sizeof(list<RAT::DBLink*>) );
      instance.SetNew(&new_listlERATcLcLDBLinkmUgR);
      instance.SetNewArray(&newArray_listlERATcLcLDBLinkmUgR);
      instance.SetDelete(&delete_listlERATcLcLDBLinkmUgR);
      instance.SetDeleteArray(&deleteArray_listlERATcLcLDBLinkmUgR);
      instance.SetDestructor(&destruct_listlERATcLcLDBLinkmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< list<RAT::DBLink*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const list<RAT::DBLink*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *listlERATcLcLDBLinkmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const list<RAT::DBLink*>*)0x0)->GetClass();
      listlERATcLcLDBLinkmUgR_TClassManip(theClass);
   return theClass;
   }

   static void listlERATcLcLDBLinkmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_listlERATcLcLDBLinkmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) list<RAT::DBLink*> : new list<RAT::DBLink*>;
   }
   static void *newArray_listlERATcLcLDBLinkmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) list<RAT::DBLink*>[nElements] : new list<RAT::DBLink*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_listlERATcLcLDBLinkmUgR(void *p) {
      delete ((list<RAT::DBLink*>*)p);
   }
   static void deleteArray_listlERATcLcLDBLinkmUgR(void *p) {
      delete [] ((list<RAT::DBLink*>*)p);
   }
   static void destruct_listlERATcLcLDBLinkmUgR(void *p) {
      typedef list<RAT::DBLink*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class list<RAT::DBLink*>

namespace ROOT {
   static TClass *dequelEpairlERATcLcLDBTableKeycOboolgRsPgR_Dictionary();
   static void dequelEpairlERATcLcLDBTableKeycOboolgRsPgR_TClassManip(TClass*);
   static void *new_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p = 0);
   static void *newArray_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(Long_t size, void *p);
   static void delete_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p);
   static void deleteArray_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p);
   static void destruct_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const deque<pair<RAT::DBTableKey,bool> >*)
   {
      deque<pair<RAT::DBTableKey,bool> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(deque<pair<RAT::DBTableKey,bool> >));
      static ::ROOT::TGenericClassInfo 
         instance("deque<pair<RAT::DBTableKey,bool> >", -2, "deque", 849,
                  typeid(deque<pair<RAT::DBTableKey,bool> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &dequelEpairlERATcLcLDBTableKeycOboolgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(deque<pair<RAT::DBTableKey,bool> >) );
      instance.SetNew(&new_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR);
      instance.SetNewArray(&newArray_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR);
      instance.SetDelete(&delete_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR);
      instance.SetDeleteArray(&deleteArray_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR);
      instance.SetDestructor(&destruct_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< deque<pair<RAT::DBTableKey,bool> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const deque<pair<RAT::DBTableKey,bool> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *dequelEpairlERATcLcLDBTableKeycOboolgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const deque<pair<RAT::DBTableKey,bool> >*)0x0)->GetClass();
      dequelEpairlERATcLcLDBTableKeycOboolgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void dequelEpairlERATcLcLDBTableKeycOboolgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) deque<pair<RAT::DBTableKey,bool> > : new deque<pair<RAT::DBTableKey,bool> >;
   }
   static void *newArray_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) deque<pair<RAT::DBTableKey,bool> >[nElements] : new deque<pair<RAT::DBTableKey,bool> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p) {
      delete ((deque<pair<RAT::DBTableKey,bool> >*)p);
   }
   static void deleteArray_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p) {
      delete [] ((deque<pair<RAT::DBTableKey,bool> >*)p);
   }
   static void destruct_dequelEpairlERATcLcLDBTableKeycOboolgRsPgR(void *p) {
      typedef deque<pair<RAT::DBTableKey,bool> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class deque<pair<RAT::DBTableKey,bool> >

namespace {
  void TriggerDictionaryInitialization_libRATDict_Impl() {
    static const char* headers[] = {
"RAT/DS/Root.hh",
"RAT/DS/EV.hh",
"RAT/DS/MC.hh",
"RAT/DS/MCParticle.hh",
"RAT/DS/MCPhoton.hh",
"RAT/DS/MCPMT.hh",
"RAT/DS/MCSummary.hh",
"RAT/DS/PMT.hh",
"RAT/DS/RunStore.hh",
"RAT/DS/Run.hh",
"RAT/DS/PosFit.hh",
"RAT/DS/PMTInfo.hh",
"RAT/DS/MCTrack.hh",
"RAT/DS/MCTrackStep.hh",
"RAT/DS/Calib.hh",
"RAT/DS/Centroid.hh",
"RAT/DS/PathFit.hh",
"RAT/DS/BonsaiFit.hh",
"RAT/DSReader.hh",
"RAT/DSWriter.hh",
"RAT/TrackNav.hh",
"RAT/TrackNode.hh",
"RAT/TrackCursor.hh",
"RAT/DB.hh",
"RAT/DBLink.hh",
"RAT/DBTextLoader.hh",
"RAT/DBJsonLoader.hh",
"RAT/HTTPDownloader.hh",
"RAT/Log.hh",
"RAT/ObjInt.hh",
"RAT/ObjDbl.hh",
0
    };
    static const char* includePaths[] = {
"/home/enrg/Software/WMUtils/local/include",
"/home/enrg/Software/WMUtils/local/include/Geant4",
"/home/enrg/Software/WMUtils/ratpac/build/include",
"/home/enrg/Software/WMUtils/ratpac/src/ds",
"/home/enrg/Software/WMUtils/ratpac/build/include",
"",
"/home/enrg/Software/WMUtils/local/include",
"/home/enrg/Software/WMUtils/ratpac/build/src/ds/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libRATDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$bits/char_traits.h")))  __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace std{template <typename > class __attribute__((annotate("$clingAutoload$bits/memoryfwd.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/MCTrackStep.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  MCTrackStep;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/MCPhoton.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  MCPhoton;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/MCParticle.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  MCParticle;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/MCTrack.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  MCTrack;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/MCPMT.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  MCPMT;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/PMT.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  PMT;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/MC.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  MC;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/Calib.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  Calib;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/EV.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  EV;}}
class __attribute__((annotate("$clingAutoload$TVector3.h")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  TVector3;
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/PMTInfo.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/RunStore.hh")))  PMTInfo;}}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/TrackNode.hh")))  __attribute__((annotate("$clingAutoload$RAT/TrackNav.hh")))  TrackNode;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DBTable.hh")))  __attribute__((annotate("$clingAutoload$RAT/DB.hh")))  DBTable;}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  Root;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/MCSummary.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  MCSummary;}}
namespace json{class __attribute__((annotate("$clingAutoload$RAT/json.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  Value;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/Log.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  Log;}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/PosFit.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  PosFit;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/Centroid.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  Centroid;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/BonsaiFit.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  BonsaiFit;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/PathFit.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/Root.hh")))  PathFit;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/Run.hh")))  __attribute__((annotate("$clingAutoload$RAT/DS/RunStore.hh")))  Run;}}
namespace RAT{namespace DS{class __attribute__((annotate("$clingAutoload$RAT/DS/RunStore.hh")))  RunStore;}}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DSReader.hh")))  DSReader;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DSWriter.hh")))  DSWriter;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/TrackCursor.hh")))  __attribute__((annotate("$clingAutoload$RAT/TrackNav.hh")))  TrackCursor;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/TrackNav.hh")))  TrackNav;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/HTTPDownloader.hh")))  __attribute__((annotate("$clingAutoload$RAT/DB.hh")))  HTTPDownloader;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DBFieldCallback.hh")))  __attribute__((annotate("$clingAutoload$RAT/DB.hh")))  DBFieldCallback;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DBLink.hh")))  __attribute__((annotate("$clingAutoload$RAT/DB.hh")))  DBLink;}
template <typename T> class __attribute__((annotate("$clingAutoload$RAT/smart_ptr.hpp")))  __attribute__((annotate("$clingAutoload$RAT/DB.hh")))  simple_ptr_nocopy;

namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DB.hh")))  DB;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DBTextLoader.hh")))  DBTextLoader;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/DBJsonLoader.hh")))  DBJsonLoader;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/ObjInt.hh")))  ObjInt;}
namespace RAT{class __attribute__((annotate("$clingAutoload$RAT/ObjDbl.hh")))  ObjDbl;}
namespace RAT{typedef simple_ptr_nocopy<RAT::DBLink> DBLinkPtr __attribute__((annotate("$clingAutoload$RAT/DB.hh"))) ;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libRATDict dictionary payload"

#ifndef G4_STORE_TRAJECTORY
  #define G4_STORE_TRAJECTORY 1
#endif
#ifndef G4VERBOSE
  #define G4VERBOSE 1
#endif
#ifndef G4UI_USE
  #define G4UI_USE 1
#endif
#ifndef G4VIS_USE
  #define G4VIS_USE 1
#endif
#ifndef G4INTY_USE_QT
  #define G4INTY_USE_QT 1
#endif
#ifndef G4UI_USE_QT
  #define G4UI_USE_QT 1
#endif
#ifndef G4VIS_USE_OPENGLQT
  #define G4VIS_USE_OPENGLQT 1
#endif
#ifndef G4VIS_USE_OPENGL
  #define G4VIS_USE_OPENGL 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "RAT/DS/Root.hh"
#include "RAT/DS/EV.hh"
#include "RAT/DS/MC.hh"
#include "RAT/DS/MCParticle.hh"
#include "RAT/DS/MCPhoton.hh"
#include "RAT/DS/MCPMT.hh"
#include "RAT/DS/MCSummary.hh"
#include "RAT/DS/PMT.hh"
#include "RAT/DS/RunStore.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/PosFit.hh"
#include "RAT/DS/PMTInfo.hh"
#include "RAT/DS/MCTrack.hh"
#include "RAT/DS/MCTrackStep.hh"
#include "RAT/DS/Calib.hh"
#include "RAT/DS/Centroid.hh"
#include "RAT/DS/PathFit.hh"
#include "RAT/DS/BonsaiFit.hh"
#include "RAT/DSReader.hh"
#include "RAT/DSWriter.hh"
#include "RAT/TrackNav.hh"
#include "RAT/TrackNode.hh"
#include "RAT/TrackCursor.hh"
#include "RAT/DB.hh"
#include "RAT/DBLink.hh"
#include "RAT/DBTextLoader.hh"
#include "RAT/DBJsonLoader.hh"
#include "RAT/HTTPDownloader.hh"
#include "RAT/Log.hh"
#include "RAT/ObjInt.hh"
#include "RAT/ObjDbl.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RAT::DB", payloadCode, "@",
"RAT::DBFieldCallback", payloadCode, "@",
"RAT::DBJsonLoader", payloadCode, "@",
"RAT::DBLink", payloadCode, "@",
"RAT::DBLinkPtr", payloadCode, "@",
"RAT::DBTable", payloadCode, "@",
"RAT::DBTextLoader", payloadCode, "@",
"RAT::DS::BonsaiFit", payloadCode, "@",
"RAT::DS::Calib", payloadCode, "@",
"RAT::DS::Centroid", payloadCode, "@",
"RAT::DS::EV", payloadCode, "@",
"RAT::DS::MC", payloadCode, "@",
"RAT::DS::MCPMT", payloadCode, "@",
"RAT::DS::MCParticle", payloadCode, "@",
"RAT::DS::MCPhoton", payloadCode, "@",
"RAT::DS::MCSummary", payloadCode, "@",
"RAT::DS::MCTrack", payloadCode, "@",
"RAT::DS::MCTrackStep", payloadCode, "@",
"RAT::DS::PMT", payloadCode, "@",
"RAT::DS::PMTInfo", payloadCode, "@",
"RAT::DS::PathFit", payloadCode, "@",
"RAT::DS::PosFit", payloadCode, "@",
"RAT::DS::Root", payloadCode, "@",
"RAT::DS::Run", payloadCode, "@",
"RAT::DS::RunStore", payloadCode, "@",
"RAT::DSReader", payloadCode, "@",
"RAT::DSWriter", payloadCode, "@",
"RAT::HTTPDownloader", payloadCode, "@",
"RAT::Log", payloadCode, "@",
"RAT::ObjDbl", payloadCode, "@",
"RAT::ObjInt", payloadCode, "@",
"RAT::TrackCursor", payloadCode, "@",
"RAT::TrackNav", payloadCode, "@",
"RAT::TrackNode", payloadCode, "@",
"json::Value", payloadCode, "@",
"simple_ptr_nocopy<RAT::DBLink>", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libRATDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libRATDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libRATDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libRATDict() {
  TriggerDictionaryInitialization_libRATDict_Impl();
}
