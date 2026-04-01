import json
import logging
import time
from openai import OpenAI

class PatientDataAnalyzer:
    """
    [LLM 분석기 - Chat Completions 버전]
    Assistants API 대신 Chat Completions API를 사용합니다.
    파일 업로드 없이 텍스트 데이터를 직접 프롬프트에 포함하여 분석합니다.
    이 방식은 gpt-5.1 등 최신/베타 모델 사용 시 호환성이 높습니다.
    """
    def __init__(self, config: dict):
        self.config = config
        self.api_key = self.config.get("OPENAI_API_KEY")
        self.gene_name = self.config.get("GENE_NAME", "Target Gene")
        # 사용자가 원하는 모델명 (예: gpt-5.1)을 그대로 사용
        self.model = self.config.get("OPENAI_MODEL_NAME", "gpt-4o") 
        
        if not self.api_key:
            raise ValueError("OPENAI_API_KEY가 설정되지 않았습니다.")
        self.client = OpenAI(api_key=self.api_key)

    def analyze_with_text_content(self, text_content_list: list) -> dict:
        """
        추출된 텍스트 리스트를 받아 LLM에게 분석을 요청합니다.
        """
        # 1. 텍스트 합치기
        # text_content_list에는 {"type": "text", "content": "..."} 형태의 딕셔너리들이 들어옴
        full_context = ""
        for item in text_content_list:
            full_context += f"\n\n{item['content']}"

        # 2. 시스템 프롬프트 구성
        system_prompt = f"""
        You are an expert biomedical data extractor.
        Target Gene: "{self.gene_name}".
        
        Task: Analyze the provided text context (extracted from papers) to find clinical data for patients with {self.gene_name} variants.
        
        ### CRITICAL RULES ###
        1. **Target Filtering**: Extract ONLY patients with variants in **{self.gene_name}**.
        2. **Phenotypes**: Extract ONLY clinical symptoms (e.g., "Seizures"). Do NOT include variant strings in phenotypes.
        3. **Variant Extraction**:
           - **SNVs**: Construct standard HGVS format (e.g., c.123A>G, p.Arg41His) if Pos/Ref/Alt are separated.
           - **CNVs**: Keep the original description exactly (e.g., "Deletion exon 1-3").
        4. **Merge** duplicate patients found across different sections.
        
        ### OUTPUT FORMAT ###
        Return ONLY a valid JSON object. Do not add markdown blocks like ```json.
        {{
          "patients": [
            {{
              "patient_id": "string",
              "sex": "string",
              "age": "string",
              "seizure_onset": "string",
              "seizure_status": "string",
              "phenotypes": "string",
              "variant_c": "string",
              "variant_p": "string",
              "source_file": "string"
            }}
          ]
        }}
        """

        # 3. API 호출 (Chat Completions)
        for attempt in range(3):
            try:
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": f"Here is the data extracted from the paper:\n{full_context}"}
                    ],
                    temperature=0.1, # 데이터 추출이므로 창의성 낮춤
                    response_format={"type": "json_object"} # JSON 모드 강제
                )
                
                raw_text = response.choices[0].message.content
                return self._extract_json_safely(raw_text)

            except Exception as e:
                logging.warning(f"API Error (Attempt {attempt+1}) with model {self.model}: {e}")
                time.sleep(2)
        
        return {"patients": []}

    def _extract_json_safely(self, text):
        """JSON 파싱"""
        try:
            clean = text.replace("```json", "").replace("```", "").strip()
            return json.loads(clean)
        except Exception as e:
            logging.warning(f"JSON Parsing Failed: {e}")
        return {"patients": []}