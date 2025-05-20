import math
from scipy.stats import norm
from typing import List, Optional
import plotly.graph_objects as go

class SampleSizeCalculator:
    """
    Класс для расчета размера выборки в A/B-тестах с различными параметрами.
    
    Поддерживает:
    - Относительный и абсолютный MDE
    - Поправку Бонферонни
    - Непропорциональное разделение трафика
    - Визуализацию зависимостей
    """
    
    def __init__(self):
        self.default_alpha = 0.05
        self.default_power = 0.8
        self.default_two_tailed = True
        self.default_bonferroni = 1
        self.default_split_ratio = 0.5

    def calculate_sample_size_relative(
        self,
        p1: float, 
        relative_mde: float, 
        alpha: Optional[float] = None, 
        power: Optional[float] = None, 
        two_tailed: Optional[bool] = None,
        bonferroni_correction: Optional[int] = None,
        split_ratio: Optional[float] = None
    ) -> int:
        """
        Расчет размера выборки для относительного изменения конверсии.
        
        Parameters:
        -----------
        p1 : float
            Базовый уровень конверсии (0 < p1 < 1)
        relative_mde : float
            Относительное изменение (-1 < relative_mde < 1)
        alpha : float, optional
            Уровень значимости (по умолчанию 0.05)
        power : float, optional
            Мощность теста (по умолчанию 0.8)
        two_tailed : bool, optional
            Двусторонний тест (по умолчанию True)
        bonferroni_correction : int, optional
            Поправка на множественные сравнения (по умолчанию 1)
        split_ratio : float, optional
            Доля трафика в тестовую группу (по умолчанию 0.5)
            
        Returns:
        --------
        int
            Размер выборки на каждую группу
        """
        # Установка значений по умолчанию
        alpha = alpha if alpha is not None else self.default_alpha
        power = power if power is not None else self.default_power
        two_tailed = two_tailed if two_tailed is not None else self.default_two_tailed
        bonferroni_correction = bonferroni_correction if bonferroni_correction is not None else self.default_bonferroni
        split_ratio = split_ratio if split_ratio is not None else self.default_split_ratio
        
        # Валидация входных данных
        self._validate_inputs(p1, relative_mde, alpha, power, bonferroni_correction, split_ratio)
        
        # Рассчитываем p2 с учетом относительного изменения
        p2 = p1 * (1 + relative_mde)
        if p2 >= 1:
            raise ValueError("Результирующая p2 не может быть >= 1")
        
        # Средняя конверсия
        p_bar = (p1 + p2) / 2
        
        # Применяем поправку Бонферонни
        adjusted_alpha = alpha / bonferroni_correction
        
        # Критические значения Z
        z_alpha = norm.ppf(1 - adjusted_alpha / 2) if two_tailed else norm.ppf(1 - adjusted_alpha)
        z_power = norm.ppf(power)
        
        # Расчет размера выборки для равных групп
        n = (
            (z_alpha * math.sqrt(2 * p_bar * (1 - p_bar)) + 
             z_power * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2)))
        ) ** 2 / ((p2 - p1) ** 2)
        
        # Корректировка для неравных групп
        if split_ratio != 0.5:
            k = (1 - split_ratio) / split_ratio
            n = n * (1 + k) / (4 * k)
        
        return math.ceil(n)

    def calculate_sample_size_absolute(
        self,
        p1: float, 
        absolute_mde: float, 
        **kwargs
    ) -> int:
        """
        Расчет размера выборки для абсолютного изменения конверсии.
        
        Parameters:
        -----------
        p1 : float
            Базовый уровень конверсии
        absolute_mde : float
            Абсолютное изменение конверсии (p2 - p1)
            
        Returns:
        --------
        int
            Размер выборки на каждую группу
        """
        relative_mde = absolute_mde / p1
        return self.calculate_sample_size_relative(p1, relative_mde, **kwargs)

    def sample_size_plot(
        self,
        p1: float,
        mde_list: List[float],
        relative: bool = True,
        **kwargs
    ) -> None:
        """
        Визуализация зависимости размера выборки от MDE.
        
        Parameters:
        -----------
        p1 : float
            Базовый уровень конверсии
        mde_list : List[float]
            Список значений MDE
        relative : bool
            Если True - относительный MDE, False - абсолютный
        """
        if relative:
            calculate_fn = self.calculate_sample_size_relative
            x_title = 'Relative MDE'
        else:
            calculate_fn = self.calculate_sample_size_absolute
            x_title = 'Absolute MDE'
        
        sample_sizes = [calculate_fn(p1, mde, **kwargs) for mde in mde_list]
        
        split_ratio = kwargs.get('split_ratio', self.default_split_ratio)
        
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=mde_list,
                y=sample_sizes,
                name='Per variant',
                line=dict(color='blue')
            )
        )
        fig.add_trace(
            go.Scatter(
                x=mde_list,
                y=[size * (1 + (1 - split_ratio)/split_ratio) for size in sample_sizes],
                name='Total sample',
                line=dict(color='red', dash='dash')
            )
        )
        
        
        alpha = kwargs.get('alpha', self.default_alpha)
        power = kwargs.get('power', self.default_power)
        
        fig.update_layout(
            title=f"Sample Size vs MDE | Baseline: {p1*100:.1f}% | Power: {power*100:.0f}% | Alpha: {alpha*100:.1f}%",
            xaxis_title=x_title,
            yaxis_title="Sample Size",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
        )
        fig.show()

    def _validate_inputs(
        self,
        p1: float,
        mde: float,
        alpha: float,
        power: float,
        bonferroni_correction: int,
        split_ratio: float
    ) -> None:
        """Валидация входных параметров."""
        if not 0 < p1 < 1:
            raise ValueError("p1 должно быть между 0 и 1")
        if not -1 < mde < 1:
            raise ValueError("MDE должно быть между -1 и 1")
        if not 0 < alpha < 1:
            raise ValueError("alpha должно быть между 0 и 1")
        if not 0 < power < 1:
            raise ValueError("power должно быть между 0 и 1")
        if bonferroni_correction < 1:
            raise ValueError("bonferroni_correction должен быть >= 1")
        if not 0 < split_ratio < 1:
            raise ValueError("split_ratio должно быть между 0 и 1")


