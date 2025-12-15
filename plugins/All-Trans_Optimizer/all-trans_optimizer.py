from PyQt6.QtWidgets import QMessageBox
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

PLUGIN_NAME = "All-Trans Optimizer"

def run(main_window):
    """
    現在の分子のアルキル鎖（非環状C-C結合）をAll-Trans配座に整形する
    """
    # Access the current molecule via main_window.current_mol
    mol = getattr(main_window, "current_mol", None)

    if not mol:
        QMessageBox.warning(main_window, PLUGIN_NAME, "No molecule loaded.")
        return

    try:
        # 3Dコンフォマーの取得 (存在しない場合は作成しない)
        if mol.GetNumConformers() == 0:
             QMessageBox.warning(main_window, PLUGIN_NAME, "Molecule has no 3D coordinates.")
             return
        
        conf = mol.GetConformer()

        # SMARTSパターン: 炭素-炭素(非環状)-炭素-炭素
        # 中央の結合(!@)が環に含まれていない4連続の炭素を検索
        # [#6]は炭素原子を表します
        patt = Chem.MolFromSmarts("[#6]-[#6]!@[#6]-[#6]")
        matches = mol.GetSubstructMatches(patt)

        count = 0
        if matches:
            # マッチしたすべてのねじれ角を180(Trans)に設定
            # 順番によっては後続の変更が前の変更に影響を与える可能性がありますが、
            # 単純な適用でも直鎖構造には効果的です。
            for match in matches:
                idx1, idx2, idx3, idx4 = match
                
                # 二面角を180度(Trans)に設定
                rdMolTransforms.SetDihedralDeg(conf, idx1, idx2, idx3, idx4, 180.0)
                count += 1
            
            # ビューの更新
            if hasattr(main_window, "draw_molecule_3d"):
                main_window.draw_molecule_3d(mol)
            elif hasattr(main_window, "update_view"):
                main_window.update_view()
            elif hasattr(main_window, "gl_widget"):
                getattr(main_window.gl_widget, "update", lambda: None)()

            # Push undo state after modification
            if hasattr(main_window, "push_undo_state"):
                main_window.push_undo_state()

            QMessageBox.information(main_window, PLUGIN_NAME, f"Applied All-Trans to {count} torsions.")
        else:
            QMessageBox.information(main_window, PLUGIN_NAME, "No alkyl chains found.")

    except Exception as e:
        QMessageBox.critical(main_window, PLUGIN_NAME, f"Error: {str(e)}")