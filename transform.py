import fitz  # PyMuPDF
from PIL import Image
import io
import os

def converter_pdf_para_imagem_redimensionada(caminho_pdf, caminho_saida, largura_desejada=600):
    """
    Converte a primeira página de um arquivo PDF para uma imagem
    e a redimensiona para a largura especificada, mantendo a proporção.
    Retorna True se a conversão for bem-sucedida, False caso contrário.
    """
    try:
        # Abre o arquivo PDF
        documento_pdf = fitz.open(caminho_pdf)
        
        # Verifica se o PDF tem pelo menos uma página
        if len(documento_pdf) == 0:
            print(f"Aviso: O arquivo '{caminho_pdf}' está vazio ou corrompido e será ignorado.")
            return False
            
        # Seleciona a primeira página (índice 0)
        primeira_pagina = documento_pdf.load_page(0)
        
        # Renderiza a página para uma imagem (pixmap) com alta resolução
        pix = primeira_pagina.get_pixmap(dpi=300)
        bytes_imagem = pix.tobytes("png")
        
        # Cria um objeto de imagem a partir dos bytes
        imagem = Image.open(io.BytesIO(bytes_imagem))
        
        # Calcula a nova altura para manter a proporção
        largura_original, altura_original = imagem.size
        proporcao = largura_original / float(largura_desejada)
        altura_calculada = int(altura_original / proporcao)
        
        # Redimensiona a imagem
        imagem_redimensionada = imagem.resize((largura_desejada, altura_calculada), Image.Resampling.LANCZOS)
        
        # Salva a imagem no caminho de saída
        imagem_redimensionada.save(caminho_saida)
        
        return True
        
    except Exception as e:
        print(f"Erro ao converter o arquivo '{caminho_pdf}': {e}")
        return False
    finally:
        if 'documento_pdf' in locals() and documento_pdf:
            documento_pdf.close()

def processar_pasta_de_pdfs(pasta_entrada, pasta_saida, largura_imagem=600):
    """
    Processa todos os arquivos PDF de uma pasta, convertendo-os para imagens PNG
    e salvando-os em uma pasta de destino.
    """
    # 1. Cria a pasta de saída se ela não existir
    if not os.path.exists(pasta_saida):
        print(f"Criando a pasta de destino: '{pasta_saida}'")
        os.makedirs(pasta_saida)
        
    print(f"Lendo arquivos da pasta: '{pasta_entrada}'")
    
    # 2. Lista todos os arquivos na pasta de entrada
    arquivos_na_pasta = os.listdir(pasta_entrada)
    
    arquivos_convertidos = 0
    # 3. Itera sobre cada arquivo
    for nome_arquivo in arquivos_na_pasta:
        # Verifica se o arquivo tem a extensão .pdf
        if nome_arquivo.lower().endswith(".pdf"):
            # Monta o caminho completo para o arquivo de entrada
            caminho_pdf_completo = os.path.join(pasta_entrada, nome_arquivo)
            
            # Define o nome do arquivo de saída (mesmo nome, extensão .png)
            nome_base_arquivo = os.path.splitext(nome_arquivo)[0]
            nome_arquivo_saida = f"{nome_base_arquivo}.png"
            
            # Monta o caminho completo para o arquivo de saída
            caminho_saida_completo = os.path.join(pasta_saida, nome_arquivo_saida)
            
            print(f"Convertendo '{nome_arquivo}' -> '{nome_arquivo_saida}'...")
            
            # Chama a função de conversão
            if converter_pdf_para_imagem_redimensionada(caminho_pdf_completo, caminho_saida_completo, largura_imagem):
                arquivos_convertidos += 1

    print("-" * 20)
    print(f"Processo concluído. {arquivos_convertidos} arquivos convertidos.")

# --- Exemplo de Uso ---
if __name__ == "__main__":
    # Nome da pasta onde estão os seus arquivos PDF
    # O script vai procurar por uma pasta com este nome no mesmo local onde ele for executado.
    pasta_pdfs = "./img/pdfs"  
    
    # Nome da pasta onde as imagens PNG serão salvas.
    pasta_pngs = "./img/pngs"

    # Verifica se a pasta de entrada existe antes de iniciar
    if not os.path.isdir(pasta_pdfs):
        print(f"Erro: A pasta de entrada '{pasta_pdfs}' não foi encontrada.")
        print("Por favor, crie esta pasta e coloque seus arquivos PDF dentro dela.")
    else:
        processar_pasta_de_pdfs(pasta_entrada=pasta_pdfs, pasta_saida=pasta_pngs)