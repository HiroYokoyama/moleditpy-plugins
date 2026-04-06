from PyQt6.QtWidgets import QMessageBox
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

__version__="2026.04.06"
__author__="HiroYokoyama"
PLUGIN_NAME = "All-Trans Optimizer"

def run_plugin(context):
    """
    現在の分子のアルキル鎖（非環状C-C結合）をAll-Trans配座に整形する
    """
    mw = context.get_main_window()
    mol = context.current_mol

    if not mol:
        QMessageBox.warning(mw, PLUGIN_NAME, "No molecule loaded.")
        return

    try:
        # 3Dコンフォマーの取得 (存在しない場合は作成しない)
        if mol.GetNumConformers() == 0:
                QMessageBox.warning(mw, PLUGIN_NAME, "Molecule has no 3D coordinates.")
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
            
            # Push updated molecule back so 3D view redraws with new coordinates
            context.current_mol = mol
            context.refresh_3d_view()

            # Push undo state via V3 API
            context.push_undo_checkpoint()

            context.show_status_message(f"Applied All-Trans to {count} torsions.")
            QMessageBox.information(mw, PLUGIN_NAME, f"Applied All-Trans to {count} torsions.")
        else:
            QMessageBox.information(mw, PLUGIN_NAME, "No alkyl chains found.")

    except Exception as e:
        QMessageBox.critical(mw, PLUGIN_NAME, f"Error: {str(e)}")


_launch_fn = None

def initialize(context):
    global _launch_fn
    _launch_fn = lambda: run_plugin(context)

def run(mw):
    if hasattr(mw, 'host'):
        mw = mw.host
    if _launch_fn:
        _launch_fn()
