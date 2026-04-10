import hashlib
import sys

def calculate_file_sha256(file_path):
    # SHA-256のハッシュオブジェクトを作成
    sha256_hash = hashlib.sha256()
    
    # 読み込むブロックサイズ（例: 64KB）
    # メモリ効率のため、ファイルを一度に読み込まず分割して処理します
    chunk_size = 65536 

    try:
        with open(file_path, "rb") as f:
            # ファイルの終わりまでブロックごとに読み込んでハッシュを更新
            while chunk := f.read(chunk_size):
                sha256_hash.update(chunk)
                
        # 16進数文字列としてハッシュ値を返す
        return sha256_hash.hexdigest()
        
    except FileNotFoundError:
        return None
    except PermissionError:
        return "Permission Error"

if __name__ == "__main__":
    # 計算したいファイルのパスを指定してください
    target_file = "sample.txt" 
    
    # コマンドライン引数がある場合はそちらを優先
    if len(sys.argv) > 1:
        target_file = sys.argv[1]

    print(f"計算対象: {target_file}")
    hash_value = calculate_file_sha256(target_file)

    if hash_value and hash_value != "Permission Error":
        print(f"SHA-256: {hash_value}")
    elif hash_value == "Permission Error":
        print("エラー: ファイルを開く権限がありません。")
    else:
        print("エラー: ファイルが見つかりませんでした。")