from typing import List, Tuple
from scipy.stats import norm
import plotly.express as px
import numpy as np


class SampleSizeCalculator:
    def __init__(
        self,
        baseline_cr: float,
        power: float = 0.8,
        alpha: float = 0.05,
        ratio: float = 0.5,
        bonferroni_corrections: int = 1
    ):
        """
        Калькулятор размера выборки для A/B теста.

        :param baseline_cr: базовая конверсия (в долях, например 0.05)
        :param power: мощность теста (по умолчанию 0.8)
        :param alpha: уровень значимости (по умолчанию 0.05)
        :param ratio: доля трафика в тестовой группе (по умолчанию 0.5)
        :param bonferroni_corrections: количество сравнений (для поправки Бонферрони)
        """
        self.bcr = baseline_cr
        self.power = power
        self.alpha = alpha
        self.ratio = ratio
        self.h_c = bonferroni_corrections

    def compute(self, mde_percent: float) -> Tuple[float, Tuple[float, float]]:
        """
        Вычисляет минимальный размер выборки и доверительный интервал.

        :param mde_percent: минимально детектируемый эффект (в процентах)
        :return: размер выборки и доверительный интервал
        """
    
        if self.h_c < 1:
            raise ValueError("Bonferroni corrections must be >= 1")
    
        mde = self.bcr * (mde_percent / 100)
        z_beta = norm.ppf(self.power)
        z_alpha = norm.ppf(1 - (self.alpha / self.h_c) / 2)
    
        p1 = self.bcr
        p2 = self.bcr + mde
    
        # Корректный расчет объединенной дисперсии
        pooled_p = (p1 + p2 * self.ratio) / (1 + self.ratio)
        pooled_var = p1 * (1 - p1) + p2 * (1 - p2) * self.ratio
    
        base_n = ((z_alpha + z_beta) ** 2 * pooled_var) / (mde ** 2)
    
        # Более точный расчет доверительного интервала
        se = np.sqrt(p1*(1-p1)/base_n + p2*(1-p2)/(base_n*self.ratio))
        z_ci = norm.ppf(1 - self.alpha/2)
        ci = (max(base_n - z_ci * se * base_n, 0), base_n + z_ci * se * base_n)
    
        return base_n, ci


    def plot_curve(self, mde_range: List[int]) -> List[float]:
        """
        Строит график зависимости размера выборки от MDE.

        :param mde_range: список значений MDE (в процентах)
        :return: список размеров выборок
        """
        sizes = []
        ci_lowers = []
        ci_uppers = []

        for mde in mde_range:
            n, (ci_low, ci_high) = self.compute(mde)
            sizes.append(n)
            ci_lowers.append(ci_low)
            ci_uppers.append(ci_high)

        fig = px.line(
            x=mde_range,
            y=sizes,
            labels={'x': 'MDE (%)', 'y': 'Размер выборки'},
            title=f'Размер выборки vs MDE (Bonferroni: {self.h_c})'
        )

        fig.add_scatter(x=mde_range, y=ci_lowers, mode='lines', name='Нижняя граница CI', line=dict(dash='dot'))
        fig.add_scatter(x=mde_range, y=ci_uppers, mode='lines', name='Верхняя граница CI', line=dict(dash='dot'))

        fig.update_traces(mode='lines+markers')
        fig.show()

        return sizes
