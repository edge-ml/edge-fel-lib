//#pragma once

#ifndef EXTRACTIONDELEGATE_H
#define EXTRACTIONDELEGATE_H

#include "ExtractionHandler.h"
#include "Extractor.h"
#include <map>
#include <string>
#include <vector>

namespace eh {
	class ExtractionHandler;
}

namespace ed {

	class ExtractionDelegate
	{
		
	public:
		//Constructor that fills the handler map
		ExtractionDelegate() {
			fillHandlerMap();
		}

		//Function pointer type
		typedef float (eh::ExtractionHandler::*handler_func) (std::string, std::vector<float>&);

		//Function pointer map
	    typedef std::map<std::string, handler_func> handler_map;
		handler_map handlers;

		//List of features that take parameters
		static std::vector<std::string> parameterFeatures;

		//Cache of calculated values
		static std::map<std::string, float> calculated;
		static bool doCache;
		static void checkAndInsert(std::string, float);

		//Extraction helpers
		float extractOne(std::string, std::vector<float>&, std::map<std::string, float>&);
		std::vector<float> extractOneVectorial(std::string, std::vector<float>&, std::map<std::string, float>&);
		std::map<std::string, float> extractSome(std::vector<std::string>&, std::vector<float>&, std::map<std::string, float>&);
		std::map<std::string, float> extractAll(std::vector <float>&, std::map<std::string, float>&);
		std::vector<co::cd> extractSpectrum(std::vector<float>&);

	private:
		void fillHandlerMap();
	};
	

}
#endif
