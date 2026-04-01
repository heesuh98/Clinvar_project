import os
import logging
import tarfile
import zipfile
import requests
import random
import glob
from io import BytesIO
from xml.etree import ElementTree as ET
from Bio import Entrez
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# 로깅 설정
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class AdvancedDataProcessor:
    """
    PMID를 기반으로 OpenAlex 및 PMC API를 통해 논문 데이터를 수집하는 클래스.
    PDF 다운로드뿐만 아니라 부록 파일(ZIP, Excel) 처리 및 압축 해제 기능.
    """

    def __init__(self, config: dict):
        self.config = config
        self.download_dir = self.config.get("DOWNLOAD_ROOT", "downloads")
        self.raw_data_root = os.path.join(self.download_dir, "raw_data")
        os.makedirs(self.raw_data_root, exist_ok=True)
        
        Entrez.email = self.config.get("ENTREZ_EMAIL", "email@example.com")
        
        self.ua_list = [
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.1 Safari/605.1.15"
        ]
        
    def _new_session(self):
        """재시도 로직이 포함된 새로운 HTTP 세션을 생성합니다."""
        s = requests.Session()
        s.headers.update({"User-Agent": random.choice(self.ua_list)})
        retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
        s.mount("https://", HTTPAdapter(max_retries=retries))
        return s

    def download_raw_data(self, pmid: str):
        """
        주어진 PMID에 해당하는 논문 데이터를 다운로드하고 압축을 해제합니다.
        OpenAlex와 PMC 소스를 순차적으로 확인합니다.

        Args:
            pmid (str): 대상 논문의 PubMed ID.
        """
        save_dir = os.path.join(self.raw_data_root, str(pmid))
        os.makedirs(save_dir, exist_ok=True)

        logging.info(f"[Download] Processing PMID {pmid}...")
        
        main_pdf_path = os.path.join(save_dir, f"{pmid}_main.pdf")
        downloaded_main = False

        # 1. OpenAlex API 시도
        oa_pdf_url = self.get_openalex_main_pdf(pmid)
        if oa_pdf_url:
            if self._download_file(oa_pdf_url, main_pdf_path):
                downloaded_main = True

        # 2. PMC API 시도 (PDF 및 부록 패키지)
        pmcid = self._get_pmcid_from_pmid(pmid)
        if pmcid:
            oa_info = self._get_oa_file_urls(pmcid)
            
            # OpenAlex에서 실패했을 경우 PMC PDF 시도
            if not downloaded_main and oa_info.get('pdf_url'):
                if self._download_file(oa_info['pdf_url'], main_pdf_path):
                    downloaded_main = True

            # 부록 패키지(tar.gz) 다운로드 및 해제
            if oa_info.get('package_url'):
                self._extract_package_flat(oa_info['package_url'], save_dir, main_pdf_path)
            elif oa_info.get('xml_url') and not downloaded_main:
                self._download_file(oa_info['xml_url'], os.path.join(save_dir, "full_text.xml"))

        # 3. 폴더 내 ZIP 파일 처리
        zip_files = glob.glob(os.path.join(save_dir, "*.zip"))
        for zf in zip_files:
            logging.info(f"  [ZIP Extraction] Found zip file: {os.path.basename(zf)}")
            self._extract_zip_flat(zf, save_dir)

    def _extract_zip_flat(self, zip_path, save_dir):
        """
        ZIP 파일 압축을 해제하되, 폴더 구조를 무시하고 파일들을 최상위 디렉토리에 저장합니다.
        불필요한 이미지 파일 등은 필터링합니다.
        """
        junk_exts = {'.jpg', '.jpeg', '.png', '.gif', '.tif', '.tiff', '.eps'}
        junk_kws = ['__macosx', 'thumbs.db'] 
        
        try:
            with zipfile.ZipFile(zip_path, 'r') as z:
                for file_info in z.infolist():
                    if file_info.is_dir(): continue
                    
                    # 파일명 인코딩 보정
                    try:
                        filename = file_info.filename.encode('cp437').decode('euc-kr')
                    except:
                        filename = file_info.filename

                    base_name = os.path.basename(filename)
                    if not base_name: continue
                    
                    lower_name = base_name.lower()
                    _, ext = os.path.splitext(lower_name)

                    if any(k in lower_name for k in junk_kws): continue
                    if ext in junk_exts: continue
                    if base_name.startswith('.'): continue

                    target_path = os.path.join(save_dir, base_name)
                    if os.path.exists(target_path):
                        continue

                    with open(target_path, "wb") as f_out:
                        f_out.write(z.read(file_info.filename))
                        
            logging.info(f"     [Success] Extracted zip: {os.path.basename(zip_path)}")
            
        except Exception as e:
            logging.warning(f"     [Error] Failed to process zip: {e}")

    def get_openalex_main_pdf(self, pmid: str):
        """OpenAlex API를 통해 논문 PDF URL을 조회합니다."""
        try:
            url = f"[https://api.openalex.org/works/pmid](https://api.openalex.org/works/pmid):{pmid}"
            r = self._new_session().get(url, timeout=10)
            if r.status_code == 200:
                data = r.json()
                best_loc = data.get('best_oa_location') or {}
                if best_loc.get('pdf_url'): return best_loc['pdf_url']
                prim_loc = data.get('primary_location') or {}
                if prim_loc.get('pdf_url'): return prim_loc['pdf_url']
        except: pass
        return None

    def _extract_package_flat(self, url, save_dir, existing_main_path):
        """PMC tar.gz 패키지를 다운로드하고 필요한 파일만 평탄화하여 압축 해제합니다."""
        junk_kws = ["license", "permission", "graphic", "logo", "disclosure"]
        main_size = os.path.getsize(existing_main_path) if os.path.exists(existing_main_path) else 0

        try:
            session = self._new_session()
            r = session.get(url, stream=True, timeout=120)
            
            with tarfile.open(fileobj=BytesIO(r.content), mode="r:gz") as tar:
                for m in tar.getmembers():
                    if not m.isfile(): continue
                    
                    fn = os.path.basename(m.name)
                    ln = fn.lower()
                    
                    if any(k in ln for k in junk_kws): continue
                    if ln.endswith(('.gif', '.jpg', '.png', '.tif', '.eps')): continue 
                    
                    # 이미 다운로드한 메인 PDF와 크기가 비슷하면 중복 다운로드 방지
                    if ln.endswith('.pdf') and main_size > 0:
                        if abs(m.size - main_size) < (main_size * 0.05): continue
                    
                    dest = os.path.join(save_dir, fn)
                    if os.path.exists(dest): continue
                        
                    with open(dest, "wb") as out:
                        out.write(tar.extractfile(m).read())
        except Exception as e:
            logging.error(f"  [Error] Failed to extract TAR package: {e}")

    def _get_pmcid_from_pmid(self, pmid: str):
        """PMID를 PMCID로 변환합니다."""
        try:
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc")
            rec = Entrez.read(handle)
            handle.close()
            if rec[0].get("LinkSetDb"): return "PMC" + rec[0]["LinkSetDb"][0]["Link"][0]["Id"]
        except: return None

    def _get_oa_file_urls(self, pmcid: str) -> dict:
        """PMC OA API를 통해 파일 URL들을 조회합니다."""
        url = f"[https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=](https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=){pmcid}"
        res = {"pdf_url": None, "package_url": None, "xml_url": None}
        try:
            r = self._new_session().get(url, timeout=15)
            root = ET.fromstring(r.content)
            for link in root.findall(".//link"):
                fmt = link.attrib.get("format")
                href = link.attrib.get("href", "").replace("ftp://", "https://")
                if fmt == "pdf": res["pdf_url"] = href
                elif fmt == "tgz": res["package_url"] = href
                elif fmt == "xml": res["xml_url"] = href
        except: pass
        return res

    def _download_file(self, url, save_path):
        """단일 파일을 다운로드합니다."""
        if os.path.exists(save_path) and os.path.getsize(save_path) > 1000: return True
        try:
            r = self._new_session().get(url, stream=True, timeout=60)
            if r.status_code == 200:
                with open(save_path, "wb") as f:
                    for c in r.iter_content(8192): f.write(c)
                return True
        except: pass
        return False